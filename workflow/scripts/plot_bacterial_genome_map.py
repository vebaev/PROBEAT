#!/usr/bin/env python3
"""Draw a multi-track bacterial genome map from Bakta and downstream results."""

import csv
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from urllib.parse import unquote

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycirclize import Circos


COLORS = {
    "CDS forward": "#4F81BD",
    "CDS reverse": "#9CC2E5",
    "ABRicate": "#D62728",
    "antiSMASH": "#9467BD",
    "BAGEL4": "#E377C2",
    "dbCAN CAZymes": "#F28E2B",
    "PGMDB probiotic genes": "#2CA02C",
    "GC content": "#38A3A5",
    "GC skew +": "#2166AC",
    "GC skew -": "#B2182B",
}


def read_fasta(path):
    records = OrderedDict()
    name = None
    chunks = []
    for line in Path(path).read_text(errors="replace").splitlines():
        if line.startswith(">"):
            if name:
                records[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name:
        records[name] = "".join(chunks).upper()
    return records


def attributes(text):
    result = {}
    for item in text.split(";"):
        key, sep, value = item.partition("=")
        if sep:
            result[key] = unquote(value)
    return result


def read_gff(path, contigs):
    cds = defaultdict(list)
    loci = {}
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[0] not in contigs or fields[2] != "CDS":
            continue
        contig_len = len(contigs[fields[0]])
        start, end = int(fields[3]) - 1, int(fields[4])
        if start < 0:
            start = 0
        if end > contig_len:
            end = contig_len
        if start >= end:
            continue
        row = (start, end, fields[6])
        cds[fields[0]].append(row)
        attrs = attributes(fields[8])
        locus = attrs.get("locus_tag") or attrs.get("ID")
        if locus:
            loci[locus] = (fields[0], start, end)
    return cds, loci


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fastani_species(ani_path, references_path):
    references = {
        row.get("assembly_accession", ""): row.get("species", "")
        for row in read_tsv(references_path)
    }
    best_ani = -1.0
    best_species = ""
    for line in Path(ani_path).read_text(errors="replace").splitlines() if Path(ani_path).exists() else []:
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            ani = float(fields[2])
        except ValueError:
            continue
        accession = Path(fields[1]).stem.removesuffix("_genomic")
        if ani > best_ani:
            best_ani = ani
            best_species = references.get(accession, "")
    if best_ani < 95:
        return ""
    # NCBI organism names may append subspecies and isolate identifiers.
    # The map label should show the confirmed binomial species plus sample ID.
    words = best_species.split()
    return " ".join(words[:2]) if len(words) >= 2 else best_species


def direct_intervals(rows, contigs, seq_key, start_key, end_key):
    found = defaultdict(list)
    for row in rows:
        raw_seqid = row.get(seq_key, "")
        seqid = raw_seqid if raw_seqid in contigs else f"contig_{raw_seqid}"
        try:
            start = int(str(row.get(start_key, "")).replace(",", "")) - 1
            end = int(str(row.get(end_key, "")).replace(",", ""))
        except ValueError:
            continue
        if seqid in contigs:
            found[seqid].append((max(0, start), min(len(contigs[seqid]), end)))
    return found


def locus_intervals(rows, loci, key, predicate=lambda row: True):
    found = defaultdict(list)
    for row in rows:
        if not predicate(row):
            continue
        locus = row.get(key, "")
        if locus in loci:
            seqid, start, end = loci[locus]
            found[seqid].append((start, end))
    return found


def antismash_intervals(directory, contigs):
    found = defaultdict(list)
    for path in sorted(Path(directory).glob("*.region*.gbk")):
        text = path.read_text(errors="replace")
        contig_match = re.search(r"^LOCUS\s+(\S+)", text, re.MULTILINE)
        start_match = re.search(r"Orig\. start\s+::\s+([\d,]+)", text)
        end_match = re.search(r"Orig\. end\s+::\s+([\d,]+)", text)
        if not (contig_match and start_match and end_match):
            continue
        raw = contig_match.group(1)
        seqid = raw if raw in contigs else f"contig_{raw}"
        if seqid in contigs:
            start = int(start_match.group(1).replace(",", ""))
            end = int(end_match.group(1).replace(",", ""))
            found[seqid].append((start, min(len(contigs[seqid]), end)))
    return found


def bagel_intervals(rows, contigs):
    found = defaultdict(list)
    for row in rows:
        raw = row.get("contig") or row.get("sequence") or row.get("AOI", "").split(".")[0]
        seqid = raw if raw in contigs else f"contig_{raw}"
        try:
            start, end = int(row.get("start", "")) - 1, int(row.get("end", ""))
        except ValueError:
            continue
        if seqid in contigs:
            found[seqid].append((max(0, start), min(len(contigs[seqid]), end)))
    return found


def gc_windows(sequence, target_windows=160):
    window = max(1000, len(sequence) // target_windows)
    xs, gc, skew = [], [], []
    for start in range(0, len(sequence), window):
        chunk = sequence[start : start + window]
        g, c = chunk.count("G"), chunk.count("C")
        xs.append(start + len(chunk) / 2)
        gc.append((g + c) / max(1, len(chunk)))
        skew.append((g - c) / max(1, g + c))
    return xs, gc, skew


def add_interval_track(sector, intervals, radius, color):
    contig_len = sector.size
    track = sector.add_track(radius, r_pad_ratio=0.10)
    track.axis(fc="#F3F4F6", ec="#D1D5DB", lw=0.25)
    for start, end in intervals:
        if start < 0:
            start = 0
        if end > contig_len:
            end = contig_len
        if end > start:
            track.rect(start, end, fc=color, ec=color, lw=0)


def main():
    if len(sys.argv) != 14:
        raise SystemExit(
            "Usage: plot_bacterial_genome_map.py SAMPLE FASTA GFF ABRICATE_CARD ABRICATE_MEGARES "
            "ANTISMASH_DIR BAGEL_TSV DBCAN_TSV PGMDB_TSV FASTANI_TSV FASTANI_REFS OUT_PNG OUT_SVG"
        )
    sample, fasta, gff, card, megares, antismash_dir, bagel, dbcan, pgmdb, ani, ani_refs, out_png, out_svg = sys.argv[1:]
    contigs = read_fasta(fasta)
    cds, loci = read_gff(gff, contigs)

    tracks = {
        "ABRicate": direct_intervals(read_tsv(card) + read_tsv(megares), contigs, "SEQUENCE", "START", "END"),
        "antiSMASH": antismash_intervals(antismash_dir, contigs),
        "BAGEL4": bagel_intervals(read_tsv(bagel), contigs),
        "dbCAN CAZymes": locus_intervals(
            read_tsv(dbcan), loci, "Gene ID", lambda row: row.get("Recommend Results", "-") not in {"", "-"}
        ),
        "PGMDB probiotic genes": locus_intervals(read_tsv(pgmdb), loci, "locus_tag"),
    }

    space = 0.8 if len(contigs) > 12 else 2
    circos = Circos({name: len(seq) for name, seq in contigs.items()}, space=space)
    largest = max(len(seq) for seq in contigs.values())
    for sector in circos.sectors:
        seqid, sequence = sector.name, contigs[sector.name]
        outer = sector.add_track((94, 100))
        outer.axis(fc="#E5E7EB", ec="#64748B", lw=0.35)
        if len(sequence) >= largest * 0.05:
            sector.text(seqid, r=105, size=6, color="#334155")
        for start, end, strand in cds.get(seqid, []):
            outer.rect(start, end, r_lim=(97, 100) if strand == "+" else (94, 97),
                       fc=COLORS["CDS forward"] if strand == "+" else COLORS["CDS reverse"], ec="none")

        for label, radius in zip(tracks, [(87, 91), (80, 84), (73, 77), (66, 70), (59, 63)]):
            add_interval_track(sector, tracks[label].get(seqid, []), radius, COLORS[label])

        xs, gc, skew = gc_windows(sequence)
        gc_track = sector.add_track((42, 55))
        gc_track.axis(fc="#F8FAFC", ec="#CBD5E1", lw=0.25)
        gc_track.line(xs, gc, vmin=0, vmax=1, color=COLORS["GC content"], lw=0.75)
        skew_track = sector.add_track((24, 38))
        skew_track.axis(fc="#F8FAFC", ec="#CBD5E1", lw=0.25)
        positive_skew = [max(0, value) for value in skew]
        negative_skew = [min(0, value) for value in skew]
        skew_track.fill_between(xs, positive_skew, 0, vmin=-1, vmax=1, color=COLORS["GC skew +"], alpha=0.75)
        skew_track.fill_between(xs, negative_skew, 0, vmin=-1, vmax=1, color=COLORS["GC skew -"], alpha=0.75)

    fig = circos.plotfig(figsize=(10, 13), dpi=220)
    # Reserve a dedicated lower band for the legend instead of overlaying it on sectors.
    fig.axes[0].set_position([0.06, 0.30, 0.88, 0.64])
    total_bp = sum(map(len, contigs.values()))
    species = fastani_species(ani, ani_refs)
    ax = fig.axes[0]
    center_x, center_y = 0.5, 0.5
    if species:
        species_words = species.split()
        ax.text(center_x, center_y + 0.035, species_words[0], transform=ax.transAxes, ha="center", va="center",
                fontsize=13, fontweight="bold", color="#0F2742")
        if len(species_words) >= 2:
            ax.text(center_x, center_y + 0.005, f"{species_words[-1]} {sample}", transform=ax.transAxes, ha="center", va="center",
                    fontsize=12, fontweight="bold", color="#0F2742")
            metrics_y = center_y - 0.035
        else:
            ax.text(center_x, center_y + 0.005, f"{species} {sample}", transform=ax.transAxes, ha="center", va="center",
                    fontsize=12, fontweight="bold", color="#0F2742")
            metrics_y = center_y - 0.035
    else:
        ax.text(center_x, center_y + 0.015, sample, transform=ax.transAxes, ha="center", va="center",
                fontsize=13, fontweight="bold", color="#0F2742")
        metrics_y = center_y - 0.030
    ax.text(center_x, metrics_y, f"{total_bp / 1_000_000:.2f} Mb | {len(contigs)} contigs",
            transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#64748B")
    handles = [Patch(facecolor=COLORS[name], edgecolor="none", label=name) for name in COLORS]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=3,
               frameon=False, fontsize=9.5, handlelength=1.8, columnspacing=2.0, labelspacing=0.7)
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, facecolor="white")
    fig.savefig(out_svg, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
