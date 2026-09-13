"""Generate JSON structures for downstream visualization."""

import json
from pathlib import Path

import pandas as pd

from .output import write_json


def build_barplot_json(df: pd.DataFrame, x_col: str, y_col: str, title: str) -> dict:
    return {
        "type": "barplot",
        "title": title,
        "x": df[x_col].tolist(),
        "y": df[y_col].tolist(),
    }


def build_heatmap_json(df: pd.DataFrame, category_col: str, count_col: str, title: str) -> dict:
    return {
        "type": "heatmap",
        "title": title,
        "categories": df[category_col].tolist(),
        "values": df[count_col].tolist(),
    }


def build_sankey_json(links: list[dict]) -> dict:
    """links: list of {source, target, value}"""
    nodes = sorted({lnk["source"] for lnk in links} | {lnk["target"] for lnk in links})
    node_idx = {n: i for i, n in enumerate(nodes)}
    return {
        "type": "sankey",
        "nodes": nodes,
        "sources": [node_idx[lnk["source"]] for lnk in links],
        "targets": [node_idx[lnk["target"]] for lnk in links],
        "values": [lnk["value"] for lnk in links],
    }


def write_all_visualizations(
    enrichment_df: pd.DataFrame,
    slim_df: pd.DataFrame,
    probiotic_df: pd.DataFrame,
    outdir: Path | str,
) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Top enriched GO barplot
    top = enrichment_df.head(20)
    if not top.empty:
        go_bar = build_barplot_json(top, "GO_ID", "Fold_Enrichment", "Top 20 Enriched GO Terms")
        write_json(go_bar, outdir / "plot_go_enrichment.json")

    # 2. GO Slim abundance barplot
    if not slim_df.empty:
        slim_bar = build_barplot_json(slim_df, "GO_Slim_Category", "Gene_Count", "GO Slim Abundance")
        write_json(slim_bar, outdir / "plot_go_slim.json")

    # 3. Probiotic heatmap
    if not probiotic_df.empty:
        heat = build_heatmap_json(probiotic_df, "GO_Slim_Category", "Gene_Count", "Probiotic Category Gene Count")
        write_json(heat, outdir / "plot_probiotic_heatmap.json")

    # 4. Sankey
    sankey_links = []
    for _, row in probiotic_df.iterrows():
        func = row["Probiotic_Function"]
        slim = row["GO_Slim_Category"]
        genes = row["Genes"].split(",")
        for g in genes:
            sankey_links.append({"source": g, "target": func, "value": 1})
            sankey_links.append({"source": func, "target": slim, "value": 1})
    if sankey_links:
        sankey = build_sankey_json(sankey_links)
        write_json(sankey, outdir / "plot_sankey.json")
