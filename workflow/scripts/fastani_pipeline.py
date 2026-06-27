#!/usr/bin/env python3
"""FastANI pipeline using RefSeq genomes from NCBI."""
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from collections import defaultdict


def get_top_species(sourmash_csv, top_n=1):
    """Extract top species name from sourmash gather results."""
    with open(sourmash_csv) as f:
        row = next(csv.DictReader(f), None)
    if not row:
        return ""
    lineage = row.get("lineage", "")
    # d__Bacteria;p__Bacillota;...;s__Lacticaseibacillus paracasei
    parts = lineage.split(";")
    species = parts[-1].strip() if parts else ""
    if species.startswith("s__"):
        species = species[3:]
    return species


def ncbi_taxon_candidates(species):
    """Return NCBI query names, relaxing GTDB placeholder species to genus."""
    candidates = [species]
    match = re.match(r"^([A-Z][A-Za-z0-9_-]+)\s+sp[.\s_-]*[A-Za-z0-9_.-]+$", species or "")
    if match:
        genus = match.group(1)
        if genus not in candidates:
            candidates.append(genus)
    return candidates


def write_unavailable_fastani(ref_dir, ready_file, reason):
    ref_dir.mkdir(parents=True, exist_ok=True)
    (ref_dir / "references.tsv").write_text(
        "species\tstrain\tlabel\tassembly_accession\tassembly_level\tchromosome_topology\tsource\n"
    )
    ready_file.write_text(f"unavailable\t0\treason={reason}\n")


def write_placeholder_heatmap(out_png, message):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.axis("off")
    ax.text(
        0.5,
        0.55,
        "FastANI unavailable",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
        color="#1d2a57",
    )
    ax.text(0.5, 0.38, message, ha="center", va="center", fontsize=10, color="#5b6472")
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_empty_fastani_table(out_file):
    Path(out_file).write_text("query\treference\tani\tfragments\ttotal_fragments\n")


def find_accessions(value):
    """Recursively extract RefSeq assembly accessions from an NCBI JSON record."""
    found = set()
    if isinstance(value, dict):
        for key, item in value.items():
            if key in {"accession", "current_accession"} and isinstance(item, str) and item.startswith("GCF_"):
                found.add(item)
            found.update(find_accessions(item))
    elif isinstance(value, list):
        for item in value:
            found.update(find_accessions(item))
    return found


def run_ncbi_command(command, timeout, attempts=4):
    """Retry transient NCBI gateway/network failures with a short backoff."""
    result = None
    for attempt in range(1, attempts + 1):
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        if result.returncode == 0:
            return result
        if attempt < attempts:
            time.sleep(attempt * 5)
    return result


def ncbi_reference_metadata(report_file):
    """Return NCBI organism and strain metadata keyed by assembly accession."""
    metadata = {}
    with open(report_file) as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            accession = record.get("accession") or record.get("currentAccession", "")
            organism = record.get("organism", {})
            species = organism.get("organismName", "")
            strain = organism.get("infraspecificNames", {}).get("strain", "")
            if not strain:
                strain = record.get("assemblyInfo", {}).get("biosample", {}).get("strain", "")
            label = f"{species} strain {strain}" if strain else species
            metadata[accession] = {
                "species": species,
                "strain": strain,
                "label": label,
                "assembly_level": record.get("assemblyInfo", {}).get("assemblyLevel", ""),
            }
    return metadata


def write_reference_manifest(out_dir, species, accessions, topology_by_accession, source_label):
    report_file = Path(out_dir) / "ncbi_package" / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    metadata = ncbi_reference_metadata(report_file)
    rows = ["species\tstrain\tlabel\tassembly_accession\tassembly_level\tchromosome_topology\tsource\n"]
    for accession in accessions:
        item = metadata.get(accession, {})
        assembly_level = item.get("assembly_level") or "unknown"
        topology = topology_by_accession.get(accession, "not_checked")
        rows.append(
            f"{item.get('species', species)}\t{item.get('strain', '')}\t"
            f"{item.get('label', species)}\t{accession}\t{assembly_level}\t{topology}\t{source_label}\n"
        )
    (Path(out_dir) / "references.tsv").write_text("".join(rows))


def has_circular_chromosome(gbff_file):
    """Return true when a GenBank assembly contains a chromosome marked circular."""
    text = Path(gbff_file).read_text(errors="replace")
    for record in text.split("\n//"):
        locus = next((line for line in record.splitlines() if line.startswith("LOCUS")), "")
        definition = " ".join(
            line.strip() for line in record.splitlines()
            if line.startswith("DEFINITION") or line.startswith("            ")
        ).lower()
        if "circular" in locus.lower() and "chromosome" in definition:
            return True
    return False


def collect_refseq_accessions(datasets_bin, species, limit, assembly_level=None):
    """Collect RefSeq accessions for a species, optionally filtered by assembly level."""
    summary_cmd = [
        datasets_bin, "summary", "genome", "taxon", species,
        "--assembly-source", "RefSeq",
        "--exclude-atypical",
        "--limit", str(limit),
        "--as-json-lines",
    ]
    if assembly_level:
        summary_cmd[5:5] = ["--assembly-level", assembly_level]

    accessions = []
    summary = None
    for attempt in range(1, 5):
        summary = run_ncbi_command(summary_cmd, timeout=300)
        if summary.returncode == 0:
            candidates = []
            for line in summary.stdout.splitlines():
                if line.strip():
                    candidates.extend(sorted(find_accessions(json.loads(line))))
            accessions = list(dict.fromkeys(candidates))[:limit]
            if accessions:
                break
        if attempt < 4:
            time.sleep(attempt * 5)
    return accessions, summary


def download_and_select_references(datasets_bin, accessions, out_dir, max_refs, require_circular):
    """Download accessions and select usable FASTA references."""
    accession_file = out_dir / "ncbi_accessions.txt"
    accession_file.write_text("\n".join(accessions) + "\n")
    archive = out_dir / "ref_genomes.zip"
    download_cmd = [
        datasets_bin, "download", "genome", "accession",
        "--inputfile", str(accession_file),
        "--include", "genome,gbff",
        "--filename", str(archive),
        "--no-progressbar",
    ]
    download = run_ncbi_command(download_cmd, timeout=900)
    if download.returncode != 0:
        raise RuntimeError(f"NCBI genome download failed: {download.stderr.strip()}")

    unpacked = out_dir / "ncbi_package"
    shutil.unpack_archive(str(archive), str(unpacked))
    refs = []
    selected_accessions = []
    topology_by_accession = {}
    for accession in accessions:
        accession_dir = unpacked / "ncbi_dataset" / "data" / accession
        candidates = sorted(accession_dir.glob("*.fna"))
        gbff_files = sorted(accession_dir.glob("*.gbff"))
        is_circular = bool(gbff_files) and any(has_circular_chromosome(path) for path in gbff_files)
        topology_by_accession[accession] = "Circular" if is_circular else "not_circular_or_unverified"
        if not candidates or (require_circular and not is_circular):
            continue
        target = out_dir / f"{accession}_genomic.fna"
        shutil.copy2(candidates[0], target)
        refs.append(target)
        selected_accessions.append(accession)
        if len(refs) == max_refs:
            break
    return refs, selected_accessions, topology_by_accession


def download_reference_genomes(species, out_dir, required_genomes=20):
    """Download RefSeq genomes from NCBI, preferring complete circular genomes."""
    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)
    datasets_bin = shutil.which("datasets")
    if not datasets_bin:
        raise RuntimeError("NCBI datasets executable is not available")

    candidate_limit = max(required_genomes * 3, 60)
    accessions, summary = collect_refseq_accessions(
        datasets_bin, species, candidate_limit, assembly_level="complete"
    )
    if accessions:
        if len(accessions) < required_genomes:
            print(
                f"NCBI returned {len(accessions)} complete RefSeq genomes for {species}; "
                f"continuing with all available complete references instead of requested {required_genomes}.",
                file=sys.stderr,
            )
        refs, selected_accessions, topology_by_accession = download_and_select_references(
            datasets_bin, accessions, out_dir, required_genomes, require_circular=True
        )
        if not refs:
            raise RuntimeError(
                f"Only {len(refs)} of {len(accessions)} downloaded Complete RefSeq assemblies "
                f"have a chromosome explicitly marked circular; at least 1 is required"
            )
        if len(refs) < required_genomes:
            print(
                f"Only {len(refs)} of {len(accessions)} downloaded Complete RefSeq assemblies "
                f"have a chromosome explicitly marked circular; continuing with available "
                f"references instead of requested {required_genomes}.",
                file=sys.stderr,
            )
        write_reference_manifest(
            out_dir, species, selected_accessions, topology_by_accession, "NCBI RefSeq complete circular"
        )
        return refs, "complete_circular"

    fallback_refs = 10
    fallback_limit = max(fallback_refs * 3, 30)
    fallback_accessions, fallback_summary = collect_refseq_accessions(
        datasets_bin, species, fallback_limit, assembly_level=None
    )
    if not fallback_accessions:
        for fallback_taxon in ncbi_taxon_candidates(species)[1:]:
            print(
                f"NCBI did not accept exact taxon {species}; trying genus-level "
                f"FastANI references from {fallback_taxon}.",
                file=sys.stderr,
            )
            accessions, summary = collect_refseq_accessions(
                datasets_bin, fallback_taxon, candidate_limit, assembly_level="complete"
            )
            if accessions:
                refs, selected_accessions, topology_by_accession = download_and_select_references(
                    datasets_bin, accessions, out_dir, required_genomes, require_circular=True
                )
                if refs:
                    write_reference_manifest(
                        out_dir,
                        fallback_taxon,
                        selected_accessions,
                        topology_by_accession,
                        f"NCBI RefSeq complete circular genus fallback from {species}",
                    )
                    return refs, "complete_circular_genus_fallback"
            fallback_accessions, fallback_summary = collect_refseq_accessions(
                datasets_bin, fallback_taxon, fallback_limit, assembly_level=None
            )
            if fallback_accessions:
                refs, selected_accessions, topology_by_accession = download_and_select_references(
                    datasets_bin, fallback_accessions, out_dir, fallback_refs, require_circular=False
                )
                if refs:
                    write_reference_manifest(
                        out_dir,
                        fallback_taxon,
                        selected_accessions,
                        topology_by_accession,
                        f"NCBI RefSeq genus fallback from {species}",
                    )
                    return refs, "refseq_genus_fallback"
        raise RuntimeError(
            f"NCBI returned 0 complete RefSeq genomes and 0 fallback RefSeq genomes for {species}. "
            f"Verify NCBI network/DNS access. Last response: "
            f"{(fallback_summary.stderr if fallback_summary else summary.stderr if summary else '').strip()}"
        )
    print(
        f"NCBI returned 0 complete RefSeq genomes for {species}; "
        f"falling back to up to {fallback_refs} available RefSeq assemblies.",
        file=sys.stderr,
    )
    refs, selected_accessions, topology_by_accession = download_and_select_references(
        datasets_bin, fallback_accessions, out_dir, fallback_refs, require_circular=False
    )
    if not refs:
        raise RuntimeError(
            f"Downloaded {len(fallback_accessions)} fallback RefSeq assemblies for {species}, "
            "but none contained a usable genome FASTA."
        )
    if len(refs) < fallback_refs:
        print(
            f"Only {len(refs)} fallback RefSeq assemblies were usable for {species}; continuing.",
            file=sys.stderr,
        )
    write_reference_manifest(out_dir, species, selected_accessions, topology_by_accession, "NCBI RefSeq fallback")
    return refs, "refseq_fallback"


def run_fastani(query, ref_dir, out_file):
    """Run FastANI comparing query genome to all references."""
    ref_files = sorted(Path(ref_dir).glob("*.fna"))
    if not ref_files:
        return
    ref_list = Path(ref_dir) / "ref_list.txt"
    ref_list.write_text("\n".join(str(f) for f in ref_files))
    # Find fastani binary
    fastani_bin = os.environ.get("FASTANI_BIN", "fastani")
    if not Path(fastani_bin).exists():
        for p in os.environ.get("PATH", "").split(":"):
            if Path(p, "fastani").exists():
                fastani_bin = str(Path(p, "fastani"))
                break
    cmd = [fastani_bin, "-q", str(query), "--rl", str(ref_list), "-o", str(out_file)]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=os.environ)
    if result.returncode != 0:
        raise RuntimeError(f"FastANI failed: {result.stderr.strip()}")
    return out_file


def plot_ani_heatmap(ani_file, out_png, sample_name, ref_dir):
    """Create ANIclustermap from FastANI output using matplotlib/seaborn."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    if not Path(ani_file).exists() or Path(ani_file).stat().st_size == 0:
        return

    labels = {}
    manifest = Path(ref_dir) / "references.tsv"
    if manifest.exists():
        with open(manifest) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                labels[row.get("assembly_accession", "")] = row.get("label", "")

    # Read FastANI output: query, ref, ANI, aligned_frags, total_frags
    data = defaultdict(dict)
    with open(ani_file) as f:
        for line in f:
            q, r, ani, *_ = line.strip().split("\t")
            q_short = sample_name
            accession = Path(r).stem.removesuffix("_genomic")
            ref_label = labels.get(accession) or accession
            data[q_short][ref_label] = float(ani) / 100
    if not data:
        return

    df = pd.DataFrame(data).T
    df = df.reindex(sorted(df.columns, key=lambda col: df.loc[sample_name, col], reverse=True), axis=1)

    n = len(df.columns)
    # Keep the chart compact so labels remain readable when scaled onto A4.
    figsize = (max(10, n * 0.80), 4.8)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(df, annot=True, fmt=".1%", cmap="RdYlBu_r", vmin=0.7, vmax=1.0,
                ax=ax, cbar_kws={"label": "ANI", "shrink": 0.55}, linewidths=0.7,
                square=True, annot_kws={"fontsize": 9.5, "fontweight": "bold"})
    ax.set_title(f"ANI between {sample_name} and reference genomes", fontsize=13, fontweight="bold")
    ax.set_ylabel("")
    ax.set_xlabel("NCBI RefSeq reference genomes", fontsize=10)
    plt.xticks(rotation=60, ha="right", fontsize=9)
    plt.yticks(rotation=0, fontsize=10)
    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


def download_main():
    sourmash_csv = sys.argv[2]
    ref_dir = Path(sys.argv[3])
    ready_file = Path(sys.argv[4])
    required_genomes = int(sys.argv[5])
    species = get_top_species(sourmash_csv)
    if not species:
        print("Top species: unavailable (no sourmash GTDB hit)")
        if ref_dir.exists():
            shutil.rmtree(ref_dir)
        write_unavailable_fastani(ref_dir, ready_file, "no_sourmash_species")
        return
    print(f"Top species: {species}")
    staging_dir = ref_dir.with_name(f"{ref_dir.name}.download_tmp")
    try:
        refs, reference_mode = download_reference_genomes(species, staging_dir, required_genomes)
    except Exception:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)
        raise
    if ref_dir.exists():
        shutil.rmtree(ref_dir)
    staging_dir.rename(ref_dir)
    ready_file.write_text(f"{species}\t{len(refs)}\trequested={required_genomes}\tmode={reference_mode}\n")
    print(f"Downloaded and validated {len(refs)} RefSeq genomes ({reference_mode})")


def fastani_main():
    assembly = sys.argv[2]
    sample = sys.argv[3]
    ref_dir = Path(sys.argv[4])
    ani_out = sys.argv[5]
    heatmap_out = sys.argv[6]
    required_genomes = int(sys.argv[7])
    refs = sorted(ref_dir.glob("GCF_*_genomic.fna"))
    if not refs:
        write_empty_fastani_table(ani_out)
        write_placeholder_heatmap(heatmap_out, "No sourmash species or NCBI references were available.")
        print("FastANI unavailable: no validated references available")
        return
    if len(refs) < required_genomes:
        print(
            f"FastANI found {len(refs)} validated references; requested {required_genomes}. "
            "Continuing with available references."
        )
    run_fastani(assembly, ref_dir, ani_out)
    if not Path(ani_out).exists() or Path(ani_out).stat().st_size == 0:
        write_empty_fastani_table(ani_out)
        write_placeholder_heatmap(
            heatmap_out,
            "No ANI comparisons were produced against the available NCBI references.",
        )
        print("FastANI unavailable: no ANI comparisons were produced against the available references")
        return
    plot_ani_heatmap(ani_out, heatmap_out, sample, ref_dir)
    if not Path(heatmap_out).exists() or Path(heatmap_out).stat().st_size == 0:
        write_placeholder_heatmap(
            heatmap_out,
            "FastANI completed but no heatmap could be drawn from the available comparisons.",
        )
    print(f"FastANI done -> {ani_out}")
    print(f"Heatmap -> {heatmap_out}")


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in {"download", "run"}:
        raise SystemExit(
            "Usage: fastani_pipeline.py download SOURMASH REF_DIR READY N | "
            "run ASSEMBLY SAMPLE REF_DIR ANI_TSV HEATMAP N"
        )
    if sys.argv[1] == "download":
        download_main()
    else:
        fastani_main()


if __name__ == "__main__":
    main()
