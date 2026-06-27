#!/usr/bin/env python3
"""Extract antiSMASH region data from regions.js for the PDF report."""
import re, json, sys
from pathlib import Path


def parse_antismash_regions(regions_js_path):
    """Return list of dicts: {idx, start, end, products, genes}"""
    with open(regions_js_path) as f:
        text = f.read()

    # Extract recordData
    match = re.search(r'var\s+recordData\s*=\s*(\[.*?\]);', text, re.DOTALL)
    products_map = {}
    gene_counts = {}
    if match:
        records = json.loads(match.group(1))
        for rec in records:
            for region in rec.get('regions', []):
                # antiSMASH uses anchors such as r5c1 in all_regions. The
                # sequence ID (for example contig_2) is not part of that key.
                key = region.get('anchor')
                if not key:
                    continue
                products_map[key] = region.get('products', [])
                gene_counts[key] = len(region.get('orfs', []))

    # Extract all_regions
    match = re.search(r'var\s+all_regions\s*=\s*(\{.*?\});\s*var\s+resultsData', text, re.DOTALL)
    if not match:
        return []

    all_regions = json.loads(match.group(1))
    results = []
    for key in all_regions.get('order', []):
        region = all_regions[key]
        products = products_map.get(key, ['unknown'])
        genes = gene_counts.get(key, 0)
        results.append({
            'idx': region['idx'],
            'start': region['start'],
            'end': region['end'],
            'products': ', '.join(products),
            'genes': genes,
            'length': region['end'] - region['start'] + 1,
        })
    return results


def main():
    js = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else js.with_name(js.stem + '_regions.tsv')
    regions = parse_antismash_regions(js)

    with open(out, 'w', newline='') as f:
        import csv
        w = csv.writer(f, delimiter='\t')
        w.writerow(['Region', 'Location', 'Length', 'Products', 'Genes'])
        for r in regions:
            w.writerow([
                str(r['idx']),
                f"{r['start']:,}-{r['end']:,}",
                f"{r['length']:,} bp",
                r['products'],
                str(r['genes']),
            ])
    print("Saved: %s (%d regions)" % (out, len(regions)))

if __name__ == '__main__':
    main()
