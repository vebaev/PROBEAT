"""GO enrichment analysis using GOATOOLS Fisher exact test."""

from collections import defaultdict
from typing import Optional

import numpy as np
import pandas as pd
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag

from .utils import fdr_correction


def run_enrichment(
    foreground_genes: set[str],
    background_genes: set[str],
    gene_to_go: dict[str, set[str]],
    go_dag: GODag,
    pval_threshold: float = 0.05,
    methods: Optional[list[str]] = None,
) -> pd.DataFrame:
    """
    Run GO enrichment (Fisher exact test) with propagated GO terms.

    Parameters
    ----------
    foreground_genes : set of gene IDs flagged by PGMDB
    background_genes : set of all gene IDs with GO annotations
    gene_to_go : mapping gene -> propagated GO terms
    go_dag : loaded GODag
    pval_threshold : alpha for BH correction
    methods : GOATOOLS methods (default: ['fdr_bh'])

    Returns
    -------
    DataFrame with GO enrichment results.
    """
    if methods is None:
        methods = ["fdr_bh"]

    # GOATOOLS needs association: gene -> set of GOs
    assoc = {g: set(gos) for g, gos in gene_to_go.items()}
    study = list(foreground_genes)
    population = list(background_genes)

    g = GOEnrichmentStudy(
        population,
        assoc,
        go_dag,
        methods=methods,
        pval=pval_threshold,
    )
    results = g.run_study(study)

    rows = []
    for r in results:
        if r.enrichment == "e":  # enriched in study vs population
            go_term = go_dag.get(r.GO)
            name = go_term.name if go_term else ""
            namespace = go_term.namespace if go_term else ""
            rows.append(
                {
                    "GO_ID": r.GO,
                    "GO_Name": name,
                    "Namespace": namespace,
                    "Foreground_Count": r.study_count,
                    "Background_Count": r.pop_count,
                    "Fold_Enrichment": (
                        (r.study_count / len(foreground_genes))
                        / (r.pop_count / len(background_genes))
                        if len(background_genes) > 0 and r.pop_count > 0
                        else np.nan
                    ),
                    "Pvalue": r.p_uncorrected,
                    "FDR": r.p_fdr_bh if hasattr(r, "p_fdr_bh") else np.nan,
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("FDR", na_position="last")
    return df
