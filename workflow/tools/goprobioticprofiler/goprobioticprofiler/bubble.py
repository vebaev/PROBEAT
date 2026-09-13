"""Build general and probiotic-focused GO bubble-plot tables."""

from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from .slim import PROBIOTIC_GO_SLIM
from .utils import fdr_correction


ROOT_GOS = {"GO:0008150", "GO:0003674", "GO:0005575"}
SIGNIFICANT_FDR = 0.05
SUPPORTIVE_FDR = 0.20
SUPPORTIVE_FOLD = 2.0
SUPPORTIVE_MIN_GENES = 2


def _enrichment_counts(hit_genes, foreground, background):
    """Return Fisher enrichment metrics for a set of genes."""
    foreground_hits = set(hit_genes) & foreground
    fg_hit = len(foreground_hits)
    bg_hit = len(hit_genes)
    foreground_total = len(foreground)
    background_total = len(background)
    odds, pvalue = fisher_exact([
        [fg_hit, foreground_total - fg_hit],
        [bg_hit - fg_hit, background_total - foreground_total - bg_hit + fg_hit],
    ], alternative="greater")
    fold = (
        (fg_hit / foreground_total) / (bg_hit / background_total)
        if fg_hit and foreground_total and bg_hit and background_total else 0.0
    )
    return foreground_hits, odds, pvalue, fold


def _assign_evidence_tier(row):
    """Classify evidence without weakening the primary FDR threshold."""
    if row["FDR"] <= SIGNIFICANT_FDR:
        return "Significant"
    if (
        row["FDR"] <= SUPPORTIVE_FDR
        and row["Fold_Enrichment"] >= SUPPORTIVE_FOLD
        and row["PGMDB_Gene_Count"] >= SUPPORTIVE_MIN_GENES
    ):
        return "Supportive"
    return "Detected"


def build_general_go_profile(gene_to_go, go_dag, top_per_namespace=8):
    """Summarize the most represented informative GO terms in the genome."""
    term_genes = defaultdict(set)
    for gene, terms in gene_to_go.items():
        for go_id in terms:
            term = go_dag.get(go_id)
            if term and go_id not in ROOT_GOS and getattr(term, "depth", 0) >= 2:
                term_genes[go_id].add(gene)

    total = max(len(gene_to_go), 1)
    rows = []
    for go_id, genes in term_genes.items():
        term = go_dag[go_id]
        rows.append({
            "GO_ID": go_id,
            "GO_Name": term.name,
            "Namespace": term.namespace,
            "Gene_Count": len(genes),
            "Genome_Coverage_Percent": 100 * len(genes) / total,
            "GO_Depth": term.depth,
            "Genes": ",".join(sorted(genes)),
        })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return (
        frame.sort_values(["Namespace", "Gene_Count", "GO_Depth"], ascending=[True, False, False])
        .groupby("Namespace", group_keys=False)
        .head(top_per_namespace)
        .reset_index(drop=True)
    )


def build_probiotic_go_profile(gene_to_go, foreground_genes, go_dag):
    """Quantify unique curated probiotic GO terms among PGMDB genes."""
    background = set(gene_to_go)
    background_total = max(len(background), 1)
    foreground = set(foreground_genes) & set(gene_to_go)
    foreground_total = max(len(foreground), 1)
    go_categories = defaultdict(set)
    for category, terms in PROBIOTIC_GO_SLIM.items():
        for go_id in terms:
            go_categories[go_id].add(category)
    rows = []

    for go_id, categories in go_categories.items():
        term = go_dag.get(go_id)
        if not term:
            continue
        genome_genes = {gene for gene, gos in gene_to_go.items() if go_id in gos}
        if not genome_genes:
            continue
        probiotic_genes, odds, pvalue, fold = _enrichment_counts(
            genome_genes, foreground, background
        )
        rows.append({
            "Probiotic_Categories": "; ".join(sorted(categories)),
            "GO_ID": go_id,
            "GO_Name": term.name,
            "Namespace": term.namespace,
            "Genome_Gene_Count": len(genome_genes),
            "PGMDB_Gene_Count": len(probiotic_genes),
            "Genome_Coverage_Percent": 100 * len(genome_genes) / background_total,
            "PGMDB_Coverage_Percent": 100 * len(probiotic_genes) / foreground_total if foreground else 0.0,
            "Fold_Enrichment": fold,
            "Odds_Ratio": odds,
            "Pvalue": pvalue,
            "Genome_Genes": ",".join(sorted(genome_genes)),
            "PGMDB_Genes": ",".join(sorted(probiotic_genes)),
        })

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["FDR"] = fdr_correction(frame["Pvalue"].tolist())
    frame["Significance"] = -np.log10(frame["FDR"].clip(lower=1e-12))
    frame["Evidence_Tier"] = frame.apply(_assign_evidence_tier, axis=1)
    return frame.sort_values(
        ["PGMDB_Gene_Count", "Fold_Enrichment", "Genome_Gene_Count"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def build_probiotic_category_profile(gene_to_go, foreground_genes, go_dag):
    """Aggregate curated GO evidence into pre-defined probiotic categories."""
    background = set(gene_to_go)
    foreground = set(foreground_genes) & background
    rows = []

    for category, curated_terms in PROBIOTIC_GO_SLIM.items():
        valid_terms = {go_id for go_id in curated_terms if go_id in go_dag}
        detected_terms = {
            go_id for go_id in valid_terms
            if any(go_id in gos for gos in gene_to_go.values())
        }
        genome_genes = {
            gene for gene, gos in gene_to_go.items()
            if gos & valid_terms
        }
        if not genome_genes:
            continue
        probiotic_genes, odds, pvalue, fold = _enrichment_counts(
            genome_genes, foreground, background
        )
        rows.append({
            "Probiotic_Category": category,
            "Curated_GO_Term_Count": len(valid_terms),
            "Detected_GO_Term_Count": len(detected_terms),
            "Category_Coverage_Percent": (
                100 * len(detected_terms) / len(valid_terms) if valid_terms else 0.0
            ),
            "Genome_Gene_Count": len(genome_genes),
            "PGMDB_Gene_Count": len(probiotic_genes),
            "Genome_Coverage_Percent": (
                100 * len(genome_genes) / len(background) if background else 0.0
            ),
            "PGMDB_Coverage_Percent": (
                100 * len(probiotic_genes) / len(foreground) if foreground else 0.0
            ),
            "Fold_Enrichment": fold,
            "Odds_Ratio": odds,
            "Pvalue": pvalue,
            "Detected_GO_Terms": ",".join(sorted(detected_terms)),
            "Genome_Genes": ",".join(sorted(genome_genes)),
            "PGMDB_Genes": ",".join(sorted(probiotic_genes)),
        })

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["FDR"] = fdr_correction(frame["Pvalue"].tolist())
    frame["Significance"] = -np.log10(frame["FDR"].clip(lower=1e-12))
    frame["Evidence_Tier"] = frame.apply(_assign_evidence_tier, axis=1)
    tier_order = pd.CategoricalDtype(
        ["Significant", "Supportive", "Detected"], ordered=True
    )
    frame["Evidence_Tier"] = frame["Evidence_Tier"].astype(tier_order)
    return frame.sort_values(
        ["Evidence_Tier", "Category_Coverage_Percent", "PGMDB_Gene_Count"],
        ascending=[True, False, False],
    ).reset_index(drop=True)
