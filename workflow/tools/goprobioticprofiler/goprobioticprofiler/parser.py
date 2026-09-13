"""Parse eggNOG-mapper, PGMDB, and OBO files."""

import csv
from pathlib import Path
from typing import Optional

import pandas as pd
from goatools.obo_parser import GODag


def parse_eggnog(annotations_path: str | Path) -> pd.DataFrame:
    """
    Parse eggNOG-mapper annotations file.

    Returns DataFrame with at least columns: query_name, GOs, Preferred_name, COG_category.
    """
    path = Path(annotations_path)
    if not path.exists():
        raise FileNotFoundError(f"eggNOG file not found: {path}")

    # eggNOG header comment lines start with '#'; last comment line before data is usually '#query_name ...'
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    # Find header line
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("#") and ("query_name" in line or "query\t" in line or "\tquery\t" in line):
            header_idx = i
            break

    if header_idx is None:
        # Try without comments (some versions)
        header_idx = 0
        for i, line in enumerate(lines):
            if not line.startswith("#"):
                header_idx = i
                break

    # Remove comment markers
    header = lines[header_idx].lstrip("#").strip().split("\t")
    data = []
    for line in lines[header_idx + 1 :]:
        if line.startswith("#") or not line.strip():
            continue
        data.append(line.strip().split("\t"))

    df = pd.DataFrame(data, columns=header)
    # Rename alternative column names
    if "query" in df.columns and "query_name" not in df.columns:
        df = df.rename(columns={"query": "query_name"})
    # Ensure required columns exist
    required = {"query_name", "GOs", "COG_category"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in eggNOG file: {missing}")
    # Preferred_name may be missing in some versions; use predicted_gene_name as fallback
    if "Preferred_name" not in df.columns:
        if "predicted_gene_name" in df.columns:
            df["Preferred_name"] = df["predicted_gene_name"]
        elif "predicted_protein_name" in df.columns:
            df["Preferred_name"] = df["predicted_protein_name"]
        else:
            df["Preferred_name"] = df["query_name"]
    return df


def parse_go_terms(go_string: str) -> list[str]:
    """Split GO term string (e.g. 'GO:0008150,GO:0003824') into list."""
    if not go_string or go_string == "-":
        return []
    return [g.strip() for g in go_string.split(",") if g.strip().startswith("GO:")]


def parse_pgmdb(pgmdb_path: str | Path) -> pd.DataFrame:
    """
    Parse PGMDB probiotic marker results.
    Expected columns: GeneID, Function, Category (or locus_tag, function for raw hits)
    """
    path = Path(pgmdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PGMDB file not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)
    # Support both summary (GeneID, Function, Category) and raw hits (locus_tag, function)
    if "locus_tag" in df.columns and "function" in df.columns:
        df = df.rename(columns={"locus_tag": "GeneID", "function": "Function"})
        if "Category" not in df.columns:
            df["Category"] = df["Function"]
    elif "gene" in df.columns and "function" in df.columns:
        df = df.rename(columns={"gene": "GeneID", "function": "Function"})
        if "Category" not in df.columns:
            df["Category"] = df["Function"]
    required = {"GeneID", "Function", "Category"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in PGMDB file: {missing}")
    return df


def load_go_dag(obo_path: str | Path) -> GODag:
    """Load GO ontology DAG."""
    path = Path(obo_path)
    if not path.exists():
        raise FileNotFoundError(f"OBO file not found: {path}")
    return GODag(str(path), optional_attrs={"relationship"})


def build_gene_to_go(
    eggnog_df: pd.DataFrame, go_dag: GODag
) -> dict[str, set[str]]:
    """
    Map each gene to its propagated GO terms.
    Returns {gene_id: {go_term1, go_term2, ...}}
    """
    gene_to_go: dict[str, set[str]] = {}
    for _, row in eggnog_df.iterrows():
        gene = str(row.get("query_name", "")).strip()
        if not gene:
            continue
        raw = str(row.get("GOs", "-")).strip()
        direct = parse_go_terms(raw)
        propagated: set[str] = set()
        for go in direct:
            if go in go_dag:
                term = go_dag[go]
                propagated.add(term.id)
                propagated.update(term.get_all_parents())
        if propagated:
            gene_to_go[gene] = propagated
    return gene_to_go


def build_gene_to_name(eggnog_df: pd.DataFrame) -> dict[str, str]:
    """Map gene query_name to Preferred_name (or query_name if empty)."""
    mapping: dict[str, str] = {}
    for _, row in eggnog_df.iterrows():
        gene = str(row.get("query_name", "")).strip()
        name = str(row.get("Preferred_name", "")).strip()
        mapping[gene] = name if name and name != "-" else gene
    return mapping
