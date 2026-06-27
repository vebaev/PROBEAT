#!/usr/bin/env python3
"""Draw gene cluster diagram from antiSMASH GBK region file using matplotlib."""
import sys, os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ANTISMASH_COLORS = {
    "biosynthetic": "#810e15",
    "biosynthetic-additional": "#f16d75",
    "transport": "#0043d3",
    "regulatory": "#0d6b0d",
    "resistance": "#c000c0",
    "other": "#737373",
}


def read_gbk_features(path):
    """Parse GBK file and return list of CDS features."""
    features = []
    region_product = ""
    with open(path) as f:
        in_features = False
        current = None
        for line in f:
            line = line.rstrip()
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                if "/product=" in line and line.strip().startswith("/product="):
                    region_product = line.split('"')[1]
                continue
            if line.startswith("     CDS") or line.startswith("     tRNA") or line.startswith("     rRNA"):
                current = {"type": "CDS" if line.startswith("     CDS") else line.strip().split()[0]}
                parts = line.split()
                if len(parts) >= 2:
                    loc = parts[1]
                    if "complement" in loc:
                        current["strand"] = -1
                        loc = loc.replace("complement(", "").replace(")", "")
                    else:
                        current["strand"] = 1
                    if ".." in loc:
                        s, e = loc.split("..")
                        current["start"] = int(s)
                        current["end"] = int(e)
                features.append(current)
                current = None
            elif line.startswith("     "):
                if current is not None:
                    line = line[5:].strip()
                    if line.startswith("/locus_tag="):
                        current["locus"] = line.split('"')[1]
                    elif line.startswith("/product="):
                        current["product"] = line.split('"')[1]
                    elif line.startswith("/gene_kind="):
                        current["kind"] = line.split('"')[1]
    return features, region_product


def draw_region(gbk_path, out_path):
    features, product = read_gbk_features(gbk_path)
    if not features:
        return

    max_end = max(f["end"] for f in features)
    region_len = max_end

    height = 2.5
    fig, ax = plt.subplots(figsize=(max(8, region_len / 1500), height))

    for ft in features:
        start = ft["start"]
        end = ft["end"]
        strand = ft.get("strand", 1)
        kind = ft.get("kind", "other")
        color = ANTISMASH_COLORS.get(kind, "#737373")

        if strand == 1:
            arrow_y = 0.7
        else:
            arrow_y = -0.7

        width = end - start
        if width < 1:
            width = max_end * 0.005  # minimum visible width

        arrow = mpatches.FancyArrow(
            start, arrow_y, width, 0, width=0.4, head_width=0.2,
            head_length=min(width * 0.3, 200),
            facecolor=color, edgecolor="#333333", linewidth=0.3
        )
        ax.add_patch(arrow)

        # Label: show locus tag if label is short enough
        label = ft.get("locus", "")
        if label and len(label) <= 15:
            mid = (start + end) / 2
            ax.text(mid, arrow_y + 0.35, label, ha="center", fontsize=5, rotation=90 if width < 500 else 0)

    ax.set_xlim(-region_len * 0.02, region_len * 1.02)
    ax.set_ylim(-1.5, 1.5)
    ax.axis("off")

    # Legend
    legend_handles = []
    seen = set()
    for ft in features:
        kind = ft.get("kind", "other")
        if kind not in seen:
            seen.add(kind)
            legend_handles.append(mpatches.Patch(color=ANTISMASH_COLORS.get(kind, "#737373"), label=kind))

    if legend_handles:
        ax.legend(handles=legend_handles, fontsize=6, loc="upper right", ncol=len(legend_handles),
                  framealpha=0.8, bbox_to_anchor=(1, 1.35))

    # Title
    name = gbk_path.stem
    ax.set_title(f"{name} — {product}  ({region_len:,} bp)", fontsize=10, fontweight="bold", pad=10)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    gbk = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else gbk.with_suffix(".png")
    draw_region(gbk, out)
    print("Saved: %s" % out)


if __name__ == "__main__":
    main()
