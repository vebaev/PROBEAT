#!/usr/bin/env python3
"""Render a compact, category-grouped overview of the bacterial workflow."""

from pathlib import Path
from html import escape
import subprocess

from graphviz import Digraph


GROUPS = {
    "READS": (
        "INPUT & READ QC",
        "FASTQC  |  FASTP  |  PORECHOP  |  RASUSA",
        "Raw reads to cleaned, normalized reads",
        "#e0f0fa",
        "#144d85",
        "#144d85",
    ),
    "ASSEMBLY": (
        "GENOME ASSEMBLY",
        "UNICYCLER  |  FLYE  |  MEDAKA",
        "Clean reads to assembled contigs",
        "#e0f0fa",
        "#144d85",
        "#144d85",
    ),
    "QUALITY": (
        "QUALITY & TAXONOMY",
        "QUAST  |  BUSCO  |  MASH\nSOURMASH GTDB  |  FASTANI",
        "Assembly quality and confirmed identity",
        "#e0f0fa",
        "#144d85",
        "#144d85",
    ),
    "ANNOTATION": (
        "GENOME ANNOTATION",
        "BAKTA FULL",
        "Proteins, genes and locus tags",
        "#fcedd6",
        "#ad570f",
        "#ad570f",
    ),
    "SAFETY": (
        "SAFETY & MOBILE ELEMENTS",
        "AMRFINDERPLUS  |  ABRICATE  |  VFDB\nPLASMIDFINDER  |  MOB-RECON  |  geNomad",
        "Resistance, virulence, plasmids and viruses",
        "#fae6e6:#fcedd6",
        "#9e2129",
        "#9e2129",
    ),
    "FUNCTION": (
        "FUNCTIONAL & METABOLIC ANALYSIS",
        "eggNOG-MAPPER  |  KOFAMSCAN  |  KEGGDECODER\nKEGGANOG  |  DBCAN V5  |  GENOMESPOT",
        "Pathways, CAZymes and growth traits",
        "#ede6f7",
        "#5c388f",
        "#5c388f",
    ),
    "PROBIOTIC": (
        "PROBIOTIC & GO TERM ANALYSIS",
        "PGMDB  |  GO TERM ANALYSIS\nGOPROBIOTICPROFILER  |  GUT-BRAIN MARKERS",
        "Probiotic genes, GO functions and gut-brain traits",
        "#e0f5e8",
        "#0d6b40",
        "#0d6b40",
    ),
    "BIOACTIVE": (
        "BIOACTIVE & METABOLIC CLUSTERS",
        "BAGEL4  |  ANTISMASH  |  GUTSMASH",
        "Bacteriocins, biosynthetic and metabolic clusters",
        "#ede6f7",
        "#5c388f",
        "#5c388f",
    ),
    "REPORT": (
        "INTEGRATION & REPORT",
        "PYCIRCLIZE  |  PDF + HTML REPORT",
        "Integrated genome map and complete sample report",
        "#1f2652",
        "#d8b65c",
        "#ffffff",
    ),
    "RESULTS": (
        "INTEGRATED RESULTS",
        "QUALITY  |  SAFETY  |  FUNCTION  |  PROBIOTIC",
        "All evidence streams combined",
        "#f8fafc",
        "#64748b",
        "#334155",
    ),
}


RANKS = [
    ["READS"],
    ["ASSEMBLY"],
    ["QUALITY", "ANNOTATION"],
    ["SAFETY", "FUNCTION"],
    ["BIOACTIVE", "PROBIOTIC"],
    ["RESULTS"],
    ["REPORT"],
]


EDGES = [
    ("READS", "ASSEMBLY"),
    ("ASSEMBLY", "QUALITY"),
    ("ASSEMBLY", "ANNOTATION"),
    ("ASSEMBLY", "SAFETY"),
    ("ASSEMBLY", "BIOACTIVE"),
    ("ASSEMBLY", "FUNCTION"),
    ("ANNOTATION", "FUNCTION"),
    ("FUNCTION", "PROBIOTIC"),
    ("QUALITY", "RESULTS"),
    ("SAFETY", "RESULTS"),
    ("BIOACTIVE", "RESULTS"),
    ("PROBIOTIC", "RESULTS"),
    ("RESULTS", "REPORT"),
]

LAYOUT_EDGES = [
    ("READS", "ASSEMBLY"),
    ("ASSEMBLY", "QUALITY"),
    ("ASSEMBLY", "ANNOTATION"),
    ("QUALITY", "SAFETY"),
    ("ANNOTATION", "FUNCTION"),
    ("SAFETY", "BIOACTIVE"),
    ("FUNCTION", "PROBIOTIC"),
    ("BIOACTIVE", "RESULTS"),
    ("PROBIOTIC", "RESULTS"),
    ("RESULTS", "REPORT"),
]


def label(title, tools, output):
    title_html = escape(title)
    tools_html = escape(tools).replace("\n", "<BR/>")
    output_html = escape(output)
    return (
        f"<<B><FONT POINT-SIZE='21'>{title_html}</FONT></B>"
        f"<BR/><BR/><B><FONT POINT-SIZE='14'>{tools_html}</FONT></B>"
        f"<BR/><BR/><FONT POINT-SIZE='11'>{output_html}</FONT>>"
    )


def manual_arrows():
    """Return hand-routed arrows matching the approved visual reference."""
    paths = [
        "M430,-1152 C430,-1120 430,-1092 430,-1062",
        "M390,-962 C360,-930 325,-904 300,-878",
        "M470,-962 C500,-930 535,-900 565,-871",
        "M306,-1025 C145,-940 100,-790 216,-687",
        "M306,-1040 C35,-915 -20,-640 45,-490",
        "M555,-1025 C770,-925 835,-800 645,-687",
        "M586,-771 C590,-742 600,-715 612,-687",
        "M645,-573 C648,-545 650,-520 652,-497",
        "M300,-764 C335,-570 360,-390 405,-299",
        "M300,-573 C330,-470 365,-365 415,-299",
        "M210,-390 C265,-350 320,-320 375,-299",
        "M652,-383 C610,-345 555,-315 500,-299",
        "M430,-199 C430,-165 430,-135 430,-108",
    ]
    path_tags = "".join(
        f'<path d="{path}" fill="none" stroke="#087d9b" stroke-width="2.1" '
        'stroke-linecap="round" marker-end="url(#workflow-arrow)"/>'
        for path in paths
    )
    return (
        '<defs><marker id="workflow-arrow" markerWidth="5.5" markerHeight="5.5" '
        'refX="5" refY="2.75" orient="auto" markerUnits="strokeWidth">'
        '<path d="M0,0 L5.5,2.75 L0,5.5 z" fill="#087d9b"/></marker></defs>'
        f'<g id="manual-workflow-arrows">{path_tags}</g>'
    )


def add_manual_arrows(svg_path):
    text = svg_path.read_text()
    insert_at = text.find("<!-- READS -->")
    if insert_at < 0:
        raise RuntimeError("Could not locate workflow nodes in rendered SVG")
    svg_path.write_text(text[:insert_at] + manual_arrows() + "\n" + text[insert_at:])


def main():
    out = Path("data/workflow_grouped")
    graph = Digraph("grouped_workflow", format="png")
    graph.attr(
        rankdir="TB",
        splines="curved",
        bgcolor="white",
        pad="0.34",
        nodesep="0.55",
        ranksep="0.68 equally",
        dpi="200",
        ratio="fill",
        size="7.2,10.54!",
        outputorder="edgesfirst",
    )
    graph.attr(
        "node",
        shape="rect",
        style="rounded,filled",
        fontname="Arial",
        margin="0.24,0.18",
        width="3.15",
        height="1.05",
        fixedsize="false",
    )
    graph.attr(
        "edge",
        color="#64748b",
        fontcolor="#475569",
        fontname="Arial",
        fontsize="9",
        arrowsize="0.72",
        penwidth="1.45",
    )

    for node, (title, tools, output, fill, stroke, text) in GROUPS.items():
        graph.node(
            node,
            label(title, tools, output),
            fillcolor=fill,
            color=stroke,
            fontcolor=text,
            penwidth="2.0",
        )

    for rank in RANKS:
        with graph.subgraph() as group:
            group.attr(rank="same")
            for node in rank:
                group.node(node)
            for left, right in zip(rank, rank[1:]):
                group.edge(left, right, style="invis", weight="40")

    for source, target in LAYOUT_EDGES:
        graph.edge(source, target, style="invis", weight="100")

    graph.format = "svg"
    graph.render(out.as_posix(), cleanup=True)
    svg_path = out.with_suffix(".svg")
    add_manual_arrows(svg_path)
    subprocess.run(
        ["mamba_envs/snakemake-conda-helper-2026-06-04/bin/rsvg-convert",
         "--width", "1439", "--output", str(out.with_suffix(".png")), str(svg_path)],
        check=True,
    )


if __name__ == "__main__":
    main()
