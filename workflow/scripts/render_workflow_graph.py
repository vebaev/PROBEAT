#!/usr/bin/env python3
from pathlib import Path

from graphviz import Digraph


COLORS = {
    "data": ("#e2e8f0", "#475569", "#0f172a"),
    "qc": ("#dbeafe", "#2563eb", "#172554"),
    "assembly": ("#dcfce7", "#16a34a", "#14532d"),
    "taxonomy": ("#cffafe", "#0891b2", "#164e63"),
    "annotation": ("#fef3c7", "#d97706", "#78350f"),
    "safety": ("#fee2e2", "#dc2626", "#7f1d1d"),
    "functional": ("#ede9fe", "#7c3aed", "#3b0764"),
    "probiotic": ("#fce7f3", "#db2777", "#831843"),
    "mobile": ("#d1fae5", "#059669", "#064e3b"),
    "report": ("#1f2652", "#d8b65c", "#ffffff"),
    "stream": ("#f8fafc", "#94a3b8", "#334155"),
}


NODES = {
    "FASTQ": ("RAW PAIRED FASTQ", "Input reads", "data"),
    "FASTQC": ("FASTQC", "Raw-read quality", "qc"),
    "FASTP": ("FASTP", "Adapter trimming + QC", "qc"),
    "RASUSA": ("RASUSA", "Subsample to 3M reads", "qc"),
    "READS_READY": ("QC READS", "Clean paired reads", "stream"),
    "UNICYCLER": ("UNICYCLER", "Genome assembly", "assembly"),
    "ASSEMBLY": ("ASSEMBLY FASTA", "Assembled genome", "stream"),
    "QUAST": ("QUAST", "Assembly statistics", "qc"),
    "MASH": ("MASH", "RefSeq identification", "taxonomy"),
    "SOURMASH": ("SOURMASH GTDB", "Taxonomy confirmation", "taxonomy"),
    "FASTANI": ("FASTANI", "ANI reference comparison", "taxonomy"),
    "FASTMLST": ("FASTMLST", "Sequence typing", "taxonomy"),
    "GTDBTK": ("GTDB-TK", "Database not installed", "taxonomy"),
    "BUSCO": ("BUSCO", "Genome completeness", "qc"),
    "BAKTA": ("BAKTA FULL", "Genome annotation", "annotation"),
    "ANNOTATION": ("BAKTA ANNOTATION", "Proteins + locus tags", "stream"),
    "AMR": ("AMRFINDERPLUS", "AMR genes", "safety"),
    "ABRICATE": ("ABRICATE", "CARD + MEGARes", "safety"),
    "VFDB": ("VFDB", "Virulence factors", "safety"),
    "PLASMID": ("PLASMIDFINDER", "Plasmid replicons", "mobile"),
    "GENOMAD": ("geNomad", "Plasmids + viruses", "mobile"),
    "MOBSUITE": ("MOB-SUITE", "Plasmid mobility typing", "mobile"),
    "ANTISMASH": ("ANTISMASH", "Biosynthetic clusters", "functional"),
    "BAGEL4": ("BAGEL4", "Bacteriocins + RiPPs", "functional"),
    "KOFAM": ("KOFAMSCAN", "KO annotations", "functional"),
    "KEGG": ("KEGGDECODER", "Pathway completeness", "functional"),
    "DBCAN": ("DBCAN V5", "CAZymes + CGC substrates", "functional"),
    "PGMDB": ("PGMDB", "Probiotic gene markers", "probiotic"),
    "EGGNOG": ("eggNOG-MAPPER", "Functional annotation", "functional"),
    "GOPROBIOTIC": ("GOPROBIOTICPROFILER", "General + probiotic GO bubbles", "probiotic"),
    "GUTBRAIN": ("GUT-BRAIN MARKERS", "Bakta + eggNOG markers", "probiotic"),
    "KEGGANOG": ("KEGGANOG", "eggNOG pathways", "functional"),
    "GENOMESPOT": ("GENOMESPOT", "Growth predictions", "functional"),
    "GUTSMASH": ("GUTSMASH", "Gut metabolic clusters", "functional"),
    "PYCIRCLIZE": ("PYCIRCLIZE", "Integrated genome map", "annotation"),
    "QC_RESULTS": ("QC + TAXONOMY RESULTS", "Reads, assembly and identity", "stream"),
    "SAFETY_RESULTS": ("SAFETY + MOBILE RESULTS", "Resistance, virulence and mobility", "stream"),
    "FUNCTION_RESULTS": ("FUNCTIONAL + PROBIOTIC RESULTS", "Metabolism and probiotic potential", "stream"),
    "REPORT": ("FINAL PDF + HTML REPORT", "Combined sample results", "report"),
}


RANKS = [
    ["FASTQ"],
    ["FASTQC", "FASTP"],
    ["RASUSA"],
    ["READS_READY"],
    ["UNICYCLER"],
    ["ASSEMBLY"],
    ["QUAST", "MASH", "BUSCO"],
    ["SOURMASH", "FASTANI", "BAKTA"],
    ["ANNOTATION"],
    ["AMR", "ABRICATE", "VFDB"],
    ["PLASMID", "GENOMAD", "MOBSUITE"],
    ["ANTISMASH", "BAGEL4", "GUTSMASH"],
    ["KOFAM", "DBCAN", "EGGNOG"],
    ["KEGG", "PGMDB", "GENOMESPOT"],
    ["GOPROBIOTIC", "GUTBRAIN", "KEGGANOG"],
    ["PYCIRCLIZE"],
    ["QC_RESULTS", "SAFETY_RESULTS", "FUNCTION_RESULTS"],
    ["REPORT"],
]


EDGES = [
    ("FASTQ", "FASTQC"),
    ("FASTQ", "FASTP"),
    ("FASTP", "RASUSA"),
    ("FASTQC", "READS_READY"),
    ("RASUSA", "READS_READY"),
    ("READS_READY", "UNICYCLER"),
    ("UNICYCLER", "ASSEMBLY"),
    ("ASSEMBLY", "QUAST"),
    ("ASSEMBLY", "MASH"),
    ("ASSEMBLY", "BUSCO"),
    ("ASSEMBLY", "BAKTA"),
    ("ASSEMBLY", "AMR"),
    ("ASSEMBLY", "ABRICATE"),
    ("ASSEMBLY", "VFDB"),
    ("ASSEMBLY", "PLASMID"),
    ("ASSEMBLY", "GENOMAD"),
    ("ASSEMBLY", "MOBSUITE"),
    ("ASSEMBLY", "GUTSMASH"),
    ("ASSEMBLY", "ANTISMASH"),
    ("MASH", "SOURMASH"),
    ("SOURMASH", "FASTANI"),
    ("BAKTA", "ANNOTATION"),
    ("ANNOTATION", "BAGEL4"),
    ("ANNOTATION", "KOFAM"),
    ("ANNOTATION", "DBCAN"),
    ("ANNOTATION", "PGMDB"),
    ("ANNOTATION", "EGGNOG"),
    ("ANNOTATION", "GENOMESPOT"),
    ("KOFAM", "KEGG"),
    ("EGGNOG", "GUTBRAIN"),
    ("EGGNOG", "GOPROBIOTIC"),
    ("PGMDB", "GOPROBIOTIC"),
    ("EGGNOG", "KEGGANOG"),
    ("ANNOTATION", "GUTBRAIN"),
    ("ANNOTATION", "PYCIRCLIZE"),
    ("ABRICATE", "PYCIRCLIZE"),
    ("ANTISMASH", "PYCIRCLIZE"),
    ("BAGEL4", "PYCIRCLIZE"),
    ("DBCAN", "PYCIRCLIZE"),
    ("PGMDB", "PYCIRCLIZE"),
    ("READS_READY", "QC_RESULTS"),
    ("QUAST", "QC_RESULTS"),
    ("SOURMASH", "QC_RESULTS"),
    ("FASTANI", "QC_RESULTS"),
    ("BUSCO", "QC_RESULTS"),
    ("BAKTA", "QC_RESULTS"),
    ("AMR", "SAFETY_RESULTS"),
    ("ABRICATE", "SAFETY_RESULTS"),
    ("VFDB", "SAFETY_RESULTS"),
    ("PLASMID", "SAFETY_RESULTS"),
    ("GENOMAD", "SAFETY_RESULTS"),
    ("MOBSUITE", "SAFETY_RESULTS"),
    ("ANTISMASH", "FUNCTION_RESULTS"),
    ("BAGEL4", "FUNCTION_RESULTS"),
    ("KEGG", "FUNCTION_RESULTS"),
    ("DBCAN", "FUNCTION_RESULTS"),
    ("PGMDB", "FUNCTION_RESULTS"),
    ("EGGNOG", "FUNCTION_RESULTS"),
    ("GOPROBIOTIC", "FUNCTION_RESULTS"),
    ("GUTBRAIN", "FUNCTION_RESULTS"),
    ("KEGGANOG", "FUNCTION_RESULTS"),
    ("GENOMESPOT", "FUNCTION_RESULTS"),
    ("GUTSMASH", "FUNCTION_RESULTS"),
    ("PYCIRCLIZE", "FUNCTION_RESULTS"),
    ("QC_RESULTS", "REPORT"),
    ("SAFETY_RESULTS", "REPORT"),
    ("FUNCTION_RESULTS", "REPORT"),
]


def html_label(name, description):
    return (
        f"<<B><FONT POINT-SIZE='20'>{name}</FONT></B>"
        f"<BR/><FONT POINT-SIZE='12'>{description}</FONT>>"
    )


def main():
    out = Path("data/workflow_vertical")
    graph = Digraph("workflow", format="png")
    graph.attr(
        rankdir="TB",
        splines="ortho",
        bgcolor="white",
        pad="0.10",
        nodesep="0.34",
        ranksep="0.38 equally",
        dpi="200",
        ratio="fill",
        size="7.2,10.21!",
        outputorder="edgesfirst",
    )
    graph.attr(
        "node",
        shape="rect",
        style="rounded,filled",
        fontname="Arial",
        margin="0.18,0.12",
        width="2.35",
        height="0.72",
        fixedsize="false",
    )
    graph.attr(
        "edge",
        color="#64748b",
        arrowsize="0.65",
        penwidth="1.15",
        tailport="s",
        headport="n",
    )

    active_nodes = {node for rank in RANKS for node in rank}
    for node, (name, description, kind) in NODES.items():
        if node not in active_nodes:
            continue
        fill, stroke, text = COLORS[kind]
        graph.node(
            node,
            html_label(name, description),
            fillcolor=fill,
            color=stroke,
            fontcolor=text,
            penwidth="1.6",
        )

    for rank in RANKS:
        with graph.subgraph() as group:
            group.attr(rank="same")
            for node in rank:
                group.node(node)
            for left, right in zip(rank, rank[1:]):
                group.edge(left, right, style="invis", weight="30")

    for upper, lower in zip(RANKS, RANKS[1:]):
        graph.edge(upper[0], lower[0], style="invis", weight="100")

    for source, target in EDGES:
        graph.edge(source, target, constraint="false")

    graph.render(out.as_posix(), cleanup=True)


if __name__ == "__main__":
    main()
