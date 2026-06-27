#!/usr/bin/env python3
"""Radar chart: COG functional categories from eggNOG annotations."""
import sys
from pathlib import Path
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COG_NAMES = {
    "A": "RNA processing",
    "B": "Chromatin",
    "C": "Energy production",
    "D": "Cell cycle",
    "E": "Amino acid transport/metabolism",
    "F": "Nucleotide transport/metabolism",
    "G": "Carbohydrate transport/metabolism",
    "H": "Coenzyme transport/metabolism",
    "I": "Lipid transport/metabolism",
    "J": "Translation",
    "K": "Transcription",
    "L": "Replication/recombination/repair",
    "M": "Cell wall/membrane biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification",
    "P": "Inorganic ion transport/metabolism",
    "Q": "Secondary metabolites",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons",
    "Z": "Cytoskeleton",
}


def main():
    eggnog = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else eggnog.with_name(eggnog.stem + "_cog_radar.png")

    cog_counter = Counter()
    with open(eggnog) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 7 and cols[6] not in ("", "-"):
                for ch in cols[6]:
                    if ch.isalpha() and ch in COG_NAMES:
                        cog_counter[ch] += 1

    total = sum(cog_counter.values())
    if total == 0:
        print("No COG annotations found")
        return

    # Sort by COG letter code, exclude S (Function unknown)
    sorted_cogs = sorted(cog_counter.items(), key=lambda x: x[0])
    sorted_cogs = [(c, v) for c, v in sorted_cogs if c != "S"]
    labels = [COG_NAMES.get(c, c) for c, _ in sorted_cogs]
    values = [v / total * 100 for _, v in sorted_cogs]
    counts = [v for _, v in sorted_cogs]

    n = len(values)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    vals = values + [values[0]]
    angles += [angles[0]]

    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    ax.fill(angles, vals, color="#059669", alpha=0.20)
    ax.plot(angles, vals, color="#059669", linewidth=2, marker="o", markersize=6)

    max_val = max(vals[:-1]) if vals else 1

    # Label each category individually with alternating radial offset
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([""] * len(angles[:-1]))
    ax.set_yticklabels([])

    for i, (angle, label) in enumerate(zip(angles[:-1], labels)):
        offset = 1.48 if i % 2 == 0 else 1.22
        ax.text(angle, max_val * offset, label,
                ha="center", va="center", fontsize=8, fontweight="bold")

    # Value labels
    for angle, val, count in zip(angles[:-1], values, counts):
        label = f"{val:.1f}% ({count})"
        ax.text(angle, val + max_val * 0.08, label, ha="center", fontsize=8, fontweight="bold", color="#064e3b")

    plt.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print("Saved: %s (%d genes, %d categories)" % (out, total, n))

if __name__ == "__main__":
    main()
