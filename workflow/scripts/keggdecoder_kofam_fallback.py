#!/usr/bin/env python3
"""Create KEGGDecoder-like functional completeness outputs from KOfamScan KOs."""

from __future__ import annotations

import html
import sys
from pathlib import Path


FUNCTIONS = {
    "Glycolysis / gluconeogenesis": {"K00844", "K01810", "K00850", "K01623", "K01803", "K00134", "K00927", "K01834", "K01689", "K00873"},
    "Pentose phosphate pathway": {"K00036", "K01057", "K07404", "K00615", "K01807", "K00616", "K01783"},
    "Pyruvate fermentation to lactate": {"K00016", "K03778", "K00101", "K00102"},
    "Acetate production": {"K00625", "K00925", "K01067"},
    "Amino acid biosynthesis": {"K00826", "K01733", "K01755", "K00820", "K01657", "K00838", "K00215", "K01915"},
    "B vitamin biosynthesis": {"K00788", "K00861", "K00946", "K00790", "K03147", "K01495", "K03635", "K00767", "K00606"},
    "Folate biosynthesis": {"K00796", "K00797", "K01737", "K01934", "K00287", "K00560"},
    "Riboflavin biosynthesis": {"K00793", "K00794", "K00795", "K11752", "K14652"},
    "Exopolysaccharide / capsule": {"K01991", "K01992", "K00975", "K00976", "K01104", "K00703", "K02847"},
    "Stress response and chaperones": {"K04077", "K04043", "K04044", "K04079", "K03671", "K03704", "K04748", "K03711"},
    "Oxidative stress response": {"K03781", "K04564", "K07217", "K03386", "K00518", "K00384"},
    "Bile and acid tolerance": {"K01533", "K01442", "K01497", "K02193", "K02194", "K02195", "K06147"},
    "Adhesion / sortase-associated proteins": {"K07273", "K07343", "K21449", "K12549"},
    "Carbohydrate transport": {"K02025", "K02026", "K02027", "K02777", "K02778", "K02779", "K02780"},
    "Peptidoglycan biosynthesis": {"K00790", "K02563", "K01929", "K01928", "K02548", "K01000", "K01921"},
    "CRISPR-associated defense": {"K19085", "K19117", "K19118", "K19119", "K07464", "K09951"},
}


def read_kos(path: Path) -> set[str]:
    kos: set[str] = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            for field in line.rstrip("\n").split("\t"):
                if field.startswith("K") and len(field) == 6 and field[1:].isdigit():
                    kos.add(field)
    return kos


def color(value: float) -> str:
    if value >= 0.75:
        return "#168a4a"
    if value >= 0.5:
        return "#8fbf3f"
    if value >= 0.25:
        return "#f0b94d"
    if value > 0:
        return "#dc6b4d"
    return "#d8dee7"


def write_tsv(sample: str, values: list[tuple[str, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("sample\t" + "\t".join(name for name, _ in values) + "\n")
        handle.write(sample + "\t" + "\t".join(f"{value:.3f}" for _, value in values) + "\n")


def write_svg(sample: str, values: list[tuple[str, float]], path: Path) -> None:
    width = 1200
    row_h = 34
    top = 74
    height = top + len(values) * row_h + 44
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='36' y='34' font-family='Arial, sans-serif' font-size='24' font-weight='700' fill='#10233d'>KEGG functional completeness: {html.escape(sample)}</text>",
        "<text x='36' y='58' font-family='Arial, sans-serif' font-size='14' fill='#64748b'>Fallback calculation from KOfamScan KO presence; values show detected KO coverage per functional group.</text>",
    ]
    bar_x, bar_w = 470, 560
    for index, (name, value) in enumerate(values):
        y = top + index * row_h
        parts.append(f"<text x='36' y='{y + 18}' font-family='Arial, sans-serif' font-size='13' fill='#1f2937'>{html.escape(name)}</text>")
        parts.append(f"<rect x='{bar_x}' y='{y}' width='{bar_w}' height='20' rx='3' fill='#eef2f7'/>")
        parts.append(f"<rect x='{bar_x}' y='{y}' width='{max(2, bar_w * value):.1f}' height='20' rx='3' fill='{color(value)}'/>")
        parts.append(f"<text x='{bar_x + bar_w + 14}' y='{y + 15}' font-family='Arial, sans-serif' font-size='12' font-weight='700' fill='#1f2937'>{value * 100:.0f}%</text>")
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> int:
    kofam = Path(sys.argv[1])
    sample = sys.argv[2]
    out_tsv = Path(sys.argv[3])
    out_svg = Path(sys.argv[4])
    kos = read_kos(kofam)
    values = []
    for name, required in FUNCTIONS.items():
        values.append((name, len(kos & required) / len(required)))
    values.sort(key=lambda item: (-item[1], item[0].lower()))
    write_tsv(sample, values, out_tsv)
    write_svg(sample, values, out_svg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
