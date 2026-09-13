"""GO Slim mapping and custom Probiotic GO Slim logic."""

from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd
from goatools.obo_parser import GODag

from .utils import safe_id


# ---------------------------------------------------------------------------
# Generic GO Slim (GO sub ontology top-level categories)
# ---------------------------------------------------------------------------
GO_SLIM_ROOTS = {
    "biological_process": "GO:0008150",
    "molecular_function": "GO:0003674",
    "cellular_component": "GO:0005575",
}

GENERIC_GO_SLIM = {
    "biological_process": {
        "GO:0008152": "metabolic_process",
        "GO:0009987": "cellular_process",
        "GO:0050896": "response_to_stimulus",
        "GO:0065007": "biological_regulation",
        "GO:0032502": "developmental_process",
        "GO:0002376": "immune_system_process",
        "GO:0044419": "interspecies_interaction",
        "GO:0044699": "single_organism_process",
    },
    "molecular_function": {
        "GO:0003824": "catalytic_activity",
        "GO:0005215": "transporter_activity",
        "GO:0005488": "binding",
        "GO:0060089": "molecular_transducer_activity",
        "GO:0038023": "signaling_receptor_activity",
        "GO:0001071": "nucleic_acid_binding_transcription_factor_activity",
        "GO:0005198": "structural_molecule_activity",
        "GO:0098772": "molecular_function_regulator",
    },
    "cellular_component": {
        "GO:0016020": "membrane",
        "GO:0005634": "nucleus",
        "GO:0005737": "cytoplasm",
        "GO:0005856": "cytoskeleton",
        "GO:0005829": "cytosol",
        "GO:0005794": "Golgi_apparatus",
        "GO:0005768": "endosome",
        "GO:0005773": "vacuole",
        "GO:0005777": "peroxisome",
        "GO:0005783": "endoplasmic_reticulum",
        "GO:0032991": "protein_complex",
        "GO:0009536": "plastid",
        "GO:0009579": "thylakoid",
        "GO:0005576": "extracellular_region",
        "GO:0031975": "envelope",
        "GO:0042579": "microbody",
    },
}


def map_to_generic_slim(go_term: str, go_dag: GODag) -> list[str]:
    """Map a GO term to its generic GO Slim category(s)."""
    if go_term not in go_dag:
        return []
    term = go_dag[go_term]
    slim = []
    for ns, roots in GENERIC_GO_SLIM.items():
        if term.namespace != ns:
            continue
        for root, label in roots.items():
            if root in go_dag and term.has_parent(root):
                slim.append(label)
        if not slim:
            # fallback to root name
            slim.append(ns)
    return slim


# ---------------------------------------------------------------------------
# Custom Probiotic GO Slim
# ---------------------------------------------------------------------------
PROBIOTIC_GO_SLIM = {
    "Gut Survival": {
        "GO:0009408": "response_to_heat",
        "GO:0034605": "cellular_response_to_heat",
        "GO:0006970": "response_to_osmotic_stress",
        "GO:0034606": "cellular_response_to_osmotic_stress",
        "GO:0006979": "response_to_oxidative_stress",
        "GO:0034599": "cellular_response_to_oxidative_stress",
        "GO:0006972": "hyperosmotic_response",
        "GO:0010447": "response_to_acidic_pH",
        "GO:0071468": "cellular_response_to_acidic_pH",
        "GO:1990451": "cellular_stress_response_to_acidic_pH",
        "GO:0015721": "bile_acid_and_bile_salt_transport",
        "GO:0015125": "bile_acid_transmembrane_transporter_activity",
    },
    "Carbohydrate Utilization": {
        "GO:0005975": "carbohydrate_metabolic_process",
        "GO:0004553": "hydrolase_activity_hydrolyzing_O_glycosyl_compounds",
        "GO:0016798": "hydrolase_activity_acting_on_glycosyl_bonds",
        "GO:0000016": "lactose_metabolic_process",
    },
    "Host Interaction": {
        "GO:0030159": "cellular_adhesion",
        "GO:0044406": "adhesion_of_symbiont_to_host",
        "GO:0044036": "cell wall organization",
        "GO:0009275": "gram_positive_bacterium-type_cell_wall",
        "GO:0044277": "cellular_response_to_bacterium",
        "GO:0048870": "cell_motility",
        "GO:0042710": "biofilm_formation",
        "GO:0009372": "quorum_sensing",
    },
    "EPS Production": {
        "GO:0000271": "polysaccharide_biosynthetic_process",
        "GO:0046353": "exopolysaccharide_biosynthetic_process",
        "GO:0045227": "capsule_polysaccharide_biosynthetic_process",
    },
    "Vitamin Production": {
        "GO:0009110": "vitamin_biosynthetic_process",
        "GO:0051188": "cofactor_biosynthetic_process",
        "GO:0006775": "fat-soluble_vitamin_biosynthetic_process",
        "GO:0009235": "cobalamin_biosynthetic_process",
    },
    "Gut-Brain Axis": {
        "GO:0006563": "L-serine_metabolic_process",
        "GO:0006565": "L-tryptophan_metabolic_process",
        "GO:0006540": "GABA_shunt",
        "GO:0006596": "polyamine_catabolic_process",
        "GO:0006598": "polyamine_biosynthetic_process",
    },
    "SCFA Potential": {
        "GO:0015949": "nucleobase-containing_small_molecule_interconversion",
        "GO:0019654": "acetate_metabolic_process",
        "GO:0006085": "acetyl-CoA_metabolic_process",
        "GO:0006089": "lactate_metabolic_process",
        "GO:0006090": "pyruvate_metabolic_process",
        "GO:0019543": "propionate_metabolic_process",
        "GO:0046353": "exopolysaccharide_biosynthetic_process",
        "GO:0005975": "carbohydrate_metabolic_process",
        "GO:0006113": "fermentation",
    },
    "Antioxidant Capacity": {
        "GO:0016209": "antioxidant_activity",
        "GO:0004791": "thioredoxin_disulfide_reductase_activity",
        "GO:0015035": "protein_disulfide_oxidoreductase_activity",
        "GO:0006749": "glutathione_metabolic_process",
    },
    "Immunomodulation": {
        "GO:0002376": "immune_system_process",
        "GO:0044267": "cellular_protein_metabolic_process",
        "GO:0006477": "protein_amino_acid_lipidation",
        "GO:0000032": "cell_wall_mannoprotein_biosynthetic_process",
    },
    "Antimicrobial Potential": {
        "GO:0042742": "defense_response_to_bacterium",
        "GO:0017000": "antibiotic_biosynthetic_process",
        "GO:0009253": "peptidoglycan_catabolic_process",
    },
    "Metal and Environmental Stress": {
        "GO:0046688": "response_to_copper_ion",
        "GO:0010038": "response_to_metal_ion",
        "GO:0010044": "response_to_metal_ion_stress",
        "GO:0006826": "iron_ion_transport",
        "GO:0000041": "transition_metal_ion_transport",
        "GO:0042221": "response_to_chemical",
    },
}


def build_probiotic_slim_map(go_dag: GODag) -> dict[str, list[str]]:
    """
    Build a mapping: GO_term -> list of probiotic slim categories.
    Uses direct mapping + propagated parents.
    """
    mapping: dict[str, list[str]] = {}
    for category, go_dict in PROBIOTIC_GO_SLIM.items():
        for go, _ in go_dict.items():
            if go not in go_dag:
                continue
            term = go_dag[go]
            all_ids = {term.id} | term.get_all_parents()
            for go_id in all_ids:
                mapping.setdefault(go_id, []).append(category)
    return mapping


def assign_genes_to_slim(
    gene_to_go: dict[str, set[str]],
    slim_map: dict[str, list[str]],
) -> dict[str, set[str]]:
    """
    Assign each gene to GO Slim categories.
    Returns {slim_category: {gene1, gene2, ...}}
    """
    slim_genes: dict[str, set[str]] = defaultdict(set)
    for gene, gos in gene_to_go.items():
        seen: set[str] = set()
        for go in gos:
            for cat in slim_map.get(go, []):
                if cat not in seen:
                    seen.add(cat)
                    slim_genes[cat].add(gene)
        # Also assign generic slim for GO terms not in probiotic slim
        # (optional: if a gene has no probiotic slim, skip)
    return dict(slim_genes)


def compute_slim_abundance(
    slim_genes: dict[str, set[str]],
) -> pd.DataFrame:
    """Return DataFrame for go_slim.tsv."""
    rows = []
    for cat, genes in sorted(slim_genes.items()):
        rows.append(
            {
                "GO_Slim_Category": cat,
                "Gene_Count": len(genes),
                "Enriched_Genes": ",".join(sorted(genes)),
            }
        )
    df = pd.DataFrame(rows)
    return df
