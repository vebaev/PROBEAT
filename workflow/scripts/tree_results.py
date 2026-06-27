#!/usr/bin/env python3
"""Generate a status-oriented results tree for one sample."""
import sys
import csv
from pathlib import Path


PROGRAMS = [
    ("FastQC", "qc/fastqc", [("{sample}_R1_fastqc.html", "Raw read R1 quality report"), ("{sample}_R2_fastqc.html", "Raw read R2 quality report")]),
    ("fastp", "qc/fastp", [("{sample}.html", "Interactive read QC report")]),
    ("Porechop", "qc/porechop", [("{sample}.json", "ONT adapter-trimming summary")]),
    ("Rasusa", "qc/rasusa", [("{sample}_PE.tsv", "Read subsampling status")]),
    ("Unicycler", "assembly", [("assembly.fasta", "Genome assembly FASTA")]),
    ("Flye", "flye", [("assembly.fasta", "ONT draft assembly")]),
    ("Medaka", "medaka", [("consensus.fasta", "Polished ONT assembly")]),
    ("QUAST", "quast", [("report.tsv", "Assembly quality statistics")]),
    ("Mash", "mash", [("top_hits.tsv", "RefSeq species identification")]),
    ("sourmash", "sourmash", [("gather.with-lineages.csv", "GTDB taxonomy classification")]),
    ("FastANI", "fastani", [("{sample}_fastani.tsv", "ANI comparisons"), ("{sample}_ani_heatmap.png", "ANI heatmap")]),
    ("BUSCO", "busco", [("busco_summary.txt", "Genome completeness assessment")]),
    ("Bakta full", "bakta", [("{sample}.gff3", "Genome annotation"), ("{sample}.faa", "Protein sequences"), ("{sample}.tsv", "Feature table"), ("{sample}.png", "Genome circular plot")]),
    ("AMRFinderPlus", "amrfinder", [("{sample}.tsv", "AMR gene scan")]),
    ("ABRicate CARD/MEGARes", "abricate", [("{sample}_card.tsv", "CARD AMR scan"), ("{sample}_megares.tsv", "MEGARes AMR scan")]),
    ("PlasmidFinder", "plasmidfinder", [("{sample}.tsv", "Plasmid replicon detection")]),
    ("VFDB", "vfdb", [("{sample}.tsv", "Virulence-factor scan")]),
    ("geNomad", "genomad/assembly_summary", [("assembly_plasmid_summary.tsv", "Plasmid classification"), ("assembly_virus_summary.tsv", "Virus/provirus classification")]),
    ("MOB-recon", "mobsuite", [("contig_report.txt", "Per-contig plasmid classification"), ("mobtyper_results.txt", "Reconstructed plasmid summary")]),
    ("eggNOG-mapper", "eggnog", [("{sample}.emapper.annotations", "Functional annotations"), ("{sample}_transport.tsv", "Transporter classification"), ("{sample}_ec.tsv", "Enzyme classes"), ("{sample}_cog_radar.png", "COG radar chart")]),
    ("KofamScan", "kofamscan", [("{sample}.tsv", "KEGG Orthology annotations")]),
    ("KEGGDecoder", "keggdecoder", [("pathway_completeness.tsv", "Pathway completeness"), ("pathway_completeness.svg", "Pathway heatmap")]),
    ("dbCAN v5", "dbcan", [("overview.tsv", "CAZyme annotations"), ("cgc/cgc_with_substrates.tsv", "CGC substrate predictions")]),
    ("PGMDB", "pgmdb", [("summary.tsv", "Probiotic marker summary"), ("hits.tsv", "Probiotic marker details"), ("pgmdb_functional_potential.png", "PGMDB gene coverage"), ("pgmdb_functional_potential.tsv", "PGMDB gene coverage data")]),
    ("GO Probiotic", "goprobioticprofiler", [("general_go_bubble.png", "Genome GO functional landscape"), ("probiotic_go_bubble.png", "Probiotic GO categories"), ("probiotic_go_categories.tsv", "Probiotic GO category evidence"), ("probiotic_go_terms.tsv", "Deduplicated probiotic GO term evidence")]),
    ("Gut-Brain markers", ".", [("gut_brain_markers.tsv", "Gut-brain axis marker table")]),
    ("KEGGaNOG", "kegganog", [("{sample}_pathways.tsv", "eggNOG pathway completeness")]),
    ("GenomeSPOT", "genomespot", [("{sample}.predictions.tsv", "Growth condition predictions")]),
    ("gutSMASH", "gutsmash", [("index.html", "Gut metabolic cluster report")]),
    ("BAGEL4", "bagel4", [("bacteriocin_clusters.tsv", "Bacteriocin and RiPP clusters")]),
    ("antiSMASH", "antismash", [("index.html", "Biosynthetic cluster report"), ("{sample}_regions.tsv", "Region summary")]),
    ("pyCirclize genome map", "genome_map", [("{sample}_genome_map.png", "Integrated bacterial genome map"), ("{sample}_genome_map.svg", "Editable vector genome map")]),
]


def sample_layout(sample):
    for path in (Path(f"data/results/{sample}/run/samples.tsv"), Path("workflow/config/samples.tsv")):
        if not path.exists():
            continue
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("sample") == sample:
                    return row.get("layout", "").upper()
    return ""


def programs_for_sample(sample):
    layout = sample_layout(sample)
    programs = []
    for program, relative_dir, files in PROGRAMS:
        if layout == "ONT" and program in {"fastp", "Unicycler"}:
            continue
        if layout not in {"ONT", "HYBRID"} and program in {"Porechop", "Flye", "Medaka"}:
            continue
        if layout == "HYBRID" and program in {"Flye", "Medaka"}:
            continue
        if program == "Rasusa":
            if layout == "HYBRID":
                files = [("{sample}_PE.tsv", "Illumina read subsampling status"), ("{sample}_ONT.tsv", "Long-read subsampling status")]
            else:
                files = [(f"{{sample}}_{layout}.tsv", "Read subsampling status")]
        programs.append((program, relative_dir, files))
    return programs


def resolved_files(results, sample, relative_dir, files):
    directory = results if relative_dir == "." else results / relative_dir
    return [
        (directory / pattern.format(sample=sample), description)
        for pattern, description in files
    ]


def main():
    sample = sys.argv[1]
    out_file = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    results = Path(f"data/results/{sample}")
    entries = []
    max_name = 0

    for program, relative_dir, files in programs_for_sample(sample):
        resolved = resolved_files(results, sample, relative_dir, files)
        successful = all(path.exists() and path.stat().st_size > 0 for path, _ in resolved)
        existing = [(path, description) for path, description in resolved if path.exists()]
        entries.append((program, successful, existing))
        for path, _ in existing:
            max_name = max(max_name, len(str(path.relative_to(results))))

    max_name += 4
    lines = [f"data/results/{sample}/"]
    for program, successful, existing in entries:
        mark = "✅" if successful else "❌"
        lines.append(f"+-- {mark} {program}/")
        if existing:
            for path, description in existing:
                relative = path.relative_to(results)
                lines.append(f"|   +-- {str(relative).ljust(max_name)} {description}")
        else:
            lines.append("|   +-- no expected output found")

    text = "\n".join(lines)
    if out_file:
        out_file.write_text(text, encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    main()
