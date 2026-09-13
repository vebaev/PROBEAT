"""Command-line interface for GOProbioticProfiler."""

import argparse
import sys
from pathlib import Path

from . import __version__
from .enrichment import run_enrichment
from .output import write_json, write_tsv
from .parser import build_gene_to_go, build_gene_to_name, parse_eggnog, parse_pgmdb, load_go_dag
from .probiotic_mapper import build_probiotic_summary, map_pgmdb_to_probiotic_categories
from .slim import assign_genes_to_slim, build_probiotic_slim_map, compute_slim_abundance
from .utils import ensure_dir
from .visualization import write_all_visualizations
from .bubble import (
    build_general_go_profile,
    build_probiotic_category_profile,
    build_probiotic_go_profile,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="goprobioticprofiler",
        description="GO enrichment and GO-Slim profiler for probiotic bacterial genomes.",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--eggnog", required=True, help="eggNOG-mapper annotations file (.emapper.annotations)")
    p.add_argument("--pgmdb", required=True, help="PGMDB probiotic marker results (.tsv)")
    p.add_argument("--obo", required=True, help="GO ontology file (go-basic.obo)")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--pval", type=float, default=0.05, help="P-value threshold for enrichment (default 0.05)")
    p.add_argument("--verbose", action="store_true", help="Print progress messages")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    outdir = ensure_dir(args.outdir)

    if args.verbose:
        print("[GOProbioticProfiler] Loading GO ontology ...", file=sys.stderr)
    go_dag = load_go_dag(args.obo)

    if args.verbose:
        print("[GOProbioticProfiler] Parsing eggNOG annotations ...", file=sys.stderr)
    eggnog_df = parse_eggnog(args.eggnog)
    gene_to_go = build_gene_to_go(eggnog_df, go_dag)
    gene_to_name = build_gene_to_name(eggnog_df)
    background_genes = set(gene_to_go.keys())

    if args.verbose:
        print("[GOProbioticProfiler] Parsing PGMDB markers ...", file=sys.stderr)
    pgmdb_df = parse_pgmdb(args.pgmdb)
    foreground_genes = set(pgmdb_df["GeneID"].tolist())
    # Keep only foreground genes that have GO annotations
    foreground_genes = foreground_genes & background_genes

    if args.verbose:
        print(f"[GOProbioticProfiler] Foreground: {len(foreground_genes)} / Background: {len(background_genes)}", file=sys.stderr)

    # 1. GO Enrichment
    if args.verbose:
        print("[GOProbioticProfiler] Running GO enrichment ...", file=sys.stderr)
    enrichment_df = run_enrichment(
        foreground_genes,
        background_genes,
        gene_to_go,
        go_dag,
        pval_threshold=args.pval,
    )
    write_tsv(enrichment_df, outdir / "go_enrichment.tsv")

    # Genome-wide GO landscape and curated probiotic GO potential.
    general_go_df = build_general_go_profile(gene_to_go, go_dag)
    write_tsv(general_go_df, outdir / "general_go_profile.tsv")
    probiotic_terms_df = build_probiotic_go_profile(gene_to_go, foreground_genes, go_dag)
    write_tsv(probiotic_terms_df, outdir / "probiotic_go_terms.tsv")
    probiotic_categories_df = build_probiotic_category_profile(
        gene_to_go, foreground_genes, go_dag
    )
    write_tsv(probiotic_categories_df, outdir / "probiotic_go_categories.tsv")

    # 2. GO Slim
    if args.verbose:
        print("[GOProbioticProfiler] Mapping GO Slim categories ...", file=sys.stderr)
    slim_map = build_probiotic_slim_map(go_dag)

    # Background (all genes) slim
    slim_genes_bg = assign_genes_to_slim(gene_to_go, slim_map)
    slim_df_bg = compute_slim_abundance(slim_genes_bg)
    slim_df_bg["FDR"] = None
    write_tsv(slim_df_bg, outdir / "go_slim_background.tsv")

    # Foreground (PGMDB) slim – only probiotic marker genes
    foreground_gene_to_go = {g: gene_to_go[g] for g in foreground_genes if g in gene_to_go}
    slim_genes_fg = assign_genes_to_slim(foreground_gene_to_go, slim_map)
    slim_df_fg = compute_slim_abundance(slim_genes_fg)
    slim_df_fg["FDR"] = None
    write_tsv(slim_df_fg, outdir / "go_slim.tsv")

    # 3. Probiotic GO Profile
    if args.verbose:
        print("[GOProbioticProfiler] Building probiotic profile ...", file=sys.stderr)
    probiotic_df = map_pgmdb_to_probiotic_categories(pgmdb_df, slim_genes_fg, gene_to_name)
    write_tsv(probiotic_df, outdir / "probiotic_go_profile.tsv")

    # 4. Summary JSON
    if args.verbose:
        print("[GOProbioticProfiler] Writing summary JSON ...", file=sys.stderr)
    summary = build_probiotic_summary(
        enrichment_df,
        slim_df_fg,
        probiotic_df,
        gene_to_go,
        slim_genes_fg,
        pgmdb_df,
    )
    write_json(summary, outdir / "probiotic_go_summary.json")

    # 5. Visualization JSONs
    if args.verbose:
        print("[GOProbioticProfiler] Writing visualization JSONs ...", file=sys.stderr)
    write_all_visualizations(enrichment_df, slim_df_fg, probiotic_df, outdir)

    if args.verbose:
        print(f"[GOProbioticProfiler] Done. Results in: {outdir}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
