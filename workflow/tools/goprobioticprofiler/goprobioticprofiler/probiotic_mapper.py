"""Map probiotic functions to GO Slim categories and genes."""

from collections import defaultdict

import pandas as pd

from .slim import PROBIOTIC_GO_SLIM


# Map lower-level probiotic functions from PGMDB to the high-level probiotic GO Slim categories.
FUNCTION_TO_CATEGORY = {
    # Gut Survival
    "acid stress": "Gut Survival",
    "bile resistance": "Gut Survival",
    "osmotic stress": "Gut Survival",
    "heat stress": "Gut Survival",
    "oxidative stress": "Gut Survival",
    "temperature stress": "Gut Survival",
    # Carbohydrate Utilization
    "carbohydrate metabolism": "Carbohydrate Utilization",
    "glycoside hydrolases": "Carbohydrate Utilization",
    "oligosaccharide utilization": "Carbohydrate Utilization",
    "prebiotic utilization": "Carbohydrate Utilization",
    # Host Interaction
    "adhesion": "Host Interaction",
    "biofilm formation": "Host Interaction",
    "cell surface proteins": "Host Interaction",
    "mucus interaction": "Host Interaction",
    # EPS Production
    "exopolysaccharide production": "EPS Production",
    "polysaccharide biosynthesis": "EPS Production",
    # Vitamin Production
    "vitamin biosynthesis": "Vitamin Production",
    "cofactor biosynthesis": "Vitamin Production",
    # Gut-Brain Axis
    "tryptophan metabolism": "Gut-Brain Axis",
    "glutamate metabolism": "Gut-Brain Axis",
    "GABA-related functions": "Gut-Brain Axis",
    "polyamine production": "Gut-Brain Axis",
    # SCFA Potential
    "SCFA production": "SCFA Potential",
    "acetate pathways": "SCFA Potential",
    "butyrate pathways": "SCFA Potential",
    "propionate pathways": "SCFA Potential",
    "fermentation": "SCFA Potential",
    # Antioxidant Capacity
    "antioxidant": "Antioxidant Capacity",
    "ROS detoxification": "Antioxidant Capacity",
    "thioredoxin system": "Antioxidant Capacity",
    "glutaredoxin system": "Antioxidant Capacity",
    # Immunomodulation
    "immunomodulation": "Immunomodulation",
    "cell wall modification": "Immunomodulation",
    "teichoic acid pathways": "Immunomodulation",
    "peptidoglycan modification": "Immunodulation",
    # Metal and Environmental Stress
    "heavy metal resistance": "Metal and Environmental Stress",
    "ion transport": "Metal and Environmental Stress",
    "environmental adaptation": "Metal and Environmental Stress",
    "intestinal persistence": "Gut Survival",
}


def normalize_function(func: str) -> str:
    return func.strip().lower()


def map_pgmdb_to_probiotic_categories(
    pgmdb_df: pd.DataFrame,
    slim_genes: dict[str, set[str]],
    gene_to_name: dict[str, str],
) -> pd.DataFrame:
    """
    Create probiotic_go_profile.tsv.
    Each row: Probiotic_Function, GO_Slim_Category, Gene_Count, Genes
    """
    rows = []
    # Group pgmdb by Category (Function in the schema above is lower-level; Category is higher-level)
    # Actually the user said columns: GeneID, Function, Category
    # We map based on Function and Category both.
    for _, row in pgmdb_df.iterrows():
        gene = row["GeneID"]
        func = row["Function"]
        cat = row["Category"]
        # Determine probiotic slim category
        slim_cat = None
        for key, target in FUNCTION_TO_CATEGORY.items():
            if key in normalize_function(func) or key in normalize_function(cat):
                slim_cat = target
                break
        if slim_cat is None:
            # Try to map from the PGMDB Category directly if it matches any of our slim categories
            for target in PROBIOTIC_GO_SLIM.keys():
                if target.lower() in normalize_function(cat):
                    slim_cat = target
                    break
        if slim_cat:
            rows.append(
                {
                    "Probiotic_Function": func,
                    "GO_Slim_Category": slim_cat,
                    "GeneID": gene,
                    "GeneName": gene_to_name.get(gene, gene),
                }
            )

    if not rows:
        return pd.DataFrame(columns=["Probiotic_Function", "GO_Slim_Category", "Gene_Count", "Genes"])

    df = pd.DataFrame(rows)
    grouped = (
        df.groupby(["Probiotic_Function", "GO_Slim_Category"])
        .apply(lambda x: pd.Series({"Gene_Count": x["GeneID"].nunique(), "Genes": ",".join(sorted(x["GeneName"].unique()))}))
        .reset_index()
    )
    return grouped


def build_probiotic_summary(
    enrichment_df: pd.DataFrame,
    slim_df: pd.DataFrame,
    probiotic_df: pd.DataFrame,
    gene_to_go: dict[str, set[str]],
    slim_genes: dict[str, set[str]],
    pgmdb_df: pd.DataFrame,
) -> dict:
    """Build JSON structure for downstream plotting."""
    top_go = (
        enrichment_df.head(20)
        .loc[:, ["GO_ID", "GO_Name", "Namespace", "Fold_Enrichment", "FDR"]]
        .to_dict(orient="records")
        if not enrichment_df.empty else []
    )
    slim_abundance = (
        slim_df.loc[:, ["GO_Slim_Category", "Gene_Count"]]
        .to_dict(orient="records")
        if not slim_df.empty else []
    )
    probiotic_profile = (
        probiotic_df.loc[:, ["Probiotic_Function", "GO_Slim_Category", "Gene_Count", "Genes"]]
        .to_dict(orient="records")
        if not probiotic_df.empty else []
    )

    # Sankey data: gene -> GO -> slim -> probiotic
    sankey_nodes: set[str] = set()
    sankey_links: list[dict] = []
    for _, row in probiotic_df.iterrows():
        func = row["Probiotic_Function"]
        slim = row["GO_Slim_Category"]
        genes = row["Genes"].split(",")
        for g in genes:
            sankey_nodes.add(g)
            sankey_nodes.add(func)
            sankey_nodes.add(slim)
            # gene -> function
            sankey_links.append({"source": g, "target": func, "value": 1})
            # function -> slim
            sankey_links.append({"source": func, "target": slim, "value": 1})

    return {
        "top_enriched_go": top_go,
        "go_slim_abundance": slim_abundance,
        "probiotic_profile": probiotic_profile,
        "sankey": {
            "nodes": list(sankey_nodes),
            "links": sankey_links,
        },
    }
