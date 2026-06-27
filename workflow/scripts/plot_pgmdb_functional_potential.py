#!/usr/bin/env python3
"""Prototype PGMDB-wide gene coverage bubble plot."""

import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


FAMILY_ORDER = [
    "Core Survival",
    "Host Interaction",
    "Metabolic Benefits",
    "Other Coverage",
]

FAMILY_COLORS = {
    "Core Survival": "#287D8E",
    "Host Interaction": "#2F855A",
    "Metabolic Benefits": "#805AD5",
    "Other Coverage": "#64748B",
}

TIER_COLORS = {
    "Strong": "#178A55",
    "Moderate": "#E69F00",
    "Detected": "#5B7FA3",
}


def split_categories(text):
    replacements = {
        "asid stress": "Acid stress",
        "temperature": "Temperature stress",
        "antioxidant": "Antioxidant",
        "bile resistance": "Bile resistance",
        "adhesions": "Adhesion",
        "heat-shock stress": "Heat stress",
    }
    categories = []
    for value in re.split(r"[/;]", str(text)):
        value = value.strip()
        if value:
            categories.append(replacements.get(value.lower(), value))
    return categories


def category_family(category):
    lowered = category.lower()
    if any(word in lowered for word in (
        "acid", "bile", "stress", "temperature", "antioxidant",
        "osmotic", "nitrite", "metal", "cold", "heat", "alkaline",
    )):
        return "Core Survival"
    if any(word in lowered for word in (
        "adhesion", "host", "persistence", "cell wall", "exopolysaccharide",
        "immunomodulation",
    )):
        return "Host Interaction"
    if any(word in lowered for word in (
        "vitamin", "scfa", "gaba", "tryptophan", "polyamine", "protease",
    )):
        return "Metabolic Benefits"
    return "Other Coverage"


def potential_tier(coverage, markers):
    if coverage >= 60 and markers >= 5:
        return "Strong"
    if coverage >= 30 and markers >= 3:
        return "Moderate"
    return "Detected"


def build_profile(database_path, hits_path):
    database = pd.read_csv(
        database_path, sep="\t", header=None,
        names=["marker", "description", "function"], dtype=str,
    ).fillna("")
    hits = pd.read_csv(hits_path, sep="\t", dtype=str).fillna("")

    database_markers = defaultdict(set)
    for row in database.itertuples(index=False):
        for category in split_categories(row.function):
            database_markers[category].add(row.marker.strip())

    hit_markers = defaultdict(set)
    hit_loci = defaultdict(set)
    bakta_markers = defaultdict(set)
    eggnog_markers = defaultdict(set)
    for row in hits.itertuples(index=False):
        for category in split_categories(row.function):
            hit_markers[category].add(row.marker.strip())
            if row.locus_tag:
                hit_loci[category].add(row.locus_tag)
            if row.source == "Bakta":
                bakta_markers[category].add(row.marker.strip())
            elif "eggNOG" in row.source:
                eggnog_markers[category].add(row.marker.strip())

    rows = []
    for category in sorted(hit_markers):
        detected = len(hit_markers[category])
        available = len(database_markers[category])
        coverage = 100 * detected / available if available else 0
        rows.append({
            "Functional_Category": category,
            "Functional_Family": category_family(category),
            "Detected_Unique_Markers": detected,
            "Database_Unique_Markers": available,
            "Database_Coverage_Percent": coverage,
            "Unique_Genome_Loci": len(hit_loci[category]),
            "Bakta_Markers": len(bakta_markers[category]),
            "eggNOG_Fallback_Markers": len(eggnog_markers[category]),
            "Potential_Tier": potential_tier(coverage, detected),
            "Detected_Markers": ",".join(sorted(hit_markers[category])),
            "Genome_Loci": ",".join(sorted(hit_loci[category])),
        })
    return pd.DataFrame(rows)


def plot_profile(profile, sample, output_path):
    profile = profile[profile["Database_Unique_Markers"] >= 3].copy()
    if profile.empty:
        return
    family_rank = {family: index for index, family in enumerate(FAMILY_ORDER)}
    profile["_family_rank"] = profile["Functional_Family"].map(family_rank)
    profile = profile.sort_values(
        ["_family_rank", "Database_Coverage_Percent", "Detected_Unique_Markers"],
        ascending=[False, True, True],
    ).reset_index(drop=True)

    y = list(range(len(profile)))
    outer_sizes = 90 + 15 * profile["Unique_Genome_Loci"].astype(float)
    inner_sizes = 25 + 7 * profile["Detected_Unique_Markers"].astype(float)
    colors = [TIER_COLORS[tier] for tier in profile["Potential_Tier"]]

    fig, ax = plt.subplots(figsize=(12, max(8, len(profile) * 0.38)))
    ax.scatter(
        profile["Database_Coverage_Percent"], y, s=outer_sizes,
        facecolors="white", edgecolors=colors, linewidths=2.2, alpha=0.95, zorder=2,
    )
    ax.scatter(
        profile["Database_Coverage_Percent"], y, s=inner_sizes,
        c=colors, edgecolors="#334155", linewidths=0.6, alpha=0.92, zorder=3,
    )

    category_labels = [
        f"{category} ({int(markers)})"
        for category, markers in zip(
            profile["Functional_Category"], profile["Database_Unique_Markers"]
        )
    ]
    ax.set_yticks(y, category_labels, fontsize=8.5)
    ax.set_xlabel("PGMDB category coverage (%)")
    ax.set_title(
        f"PGMDB Probiotic Gene Coverage: {sample}",
        fontsize=14, fontweight="bold",
    )
    ax.grid(axis="x", color="#D9E2EC", linewidth=0.7, zorder=0)
    ax.set_xlim(0, max(105, profile["Database_Coverage_Percent"].max() * 1.08))

    for boundary in range(1, len(profile)):
        if profile.loc[boundary, "Functional_Family"] != profile.loc[boundary - 1, "Functional_Family"]:
            ax.axhline(boundary - 0.5, color="#CBD5E1", linewidth=0.8)

    for family in FAMILY_ORDER:
        positions = profile.index[profile["Functional_Family"] == family].tolist()
        if positions:
            ax.text(
                1.01, sum(positions) / len(positions), family,
                transform=ax.get_yaxis_transform(), va="center", ha="left",
                fontsize=8, fontweight="bold", color=FAMILY_COLORS[family],
            )

    tier_handles = [
        plt.Line2D([], [], marker="o", linestyle="", markersize=8,
                   markerfacecolor=color, markeredgecolor="#334155", label=tier)
        for tier, color in TIER_COLORS.items()
    ]
    size_handles = [
        plt.Line2D([], [], marker="o", linestyle="", markersize=10,
                   markerfacecolor="white", markeredgecolor="#334155",
                   markeredgewidth=1.8, label="Unique genome loci (outer ring)"),
        plt.Line2D([], [], marker="o", linestyle="", markersize=6,
                   markerfacecolor="#5B7FA3", markeredgecolor="#334155",
                   label="Unique PGMDB markers (inner bubble)"),
    ]
    legend1 = ax.legend(handles=tier_handles, loc="lower right", frameon=False, fontsize=8)
    ax.add_artist(legend1)
    ax.legend(handles=size_handles, loc="upper right", frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.31, right=0.80, top=0.93, bottom=0.08)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) != 5:
        raise SystemExit(
            "Usage: plot_pgmdb_functional_potential.py DATABASE HITS OUTPUT_TSV OUTPUT_PNG"
        )
    database, hits, output_tsv, output_png = map(Path, sys.argv[1:])
    profile = build_profile(database, hits)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    profile.to_csv(output_tsv, sep="\t", index=False)
    plot_profile(profile, hits.parent.parent.name, output_png)


if __name__ == "__main__":
    main()
