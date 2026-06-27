#!/usr/bin/env python3
import csv, sys

path = sys.argv[1]
with open(path, 'r') as f:
    rdr = list(csv.reader(f, delimiter='\t'))
if not rdr:
    raise SystemExit("Empty KEGGaNOG pathway table")
hdrs = rdr[0]
row = rdr[1] if len(rdr) > 1 else ["Sample"]
keep = [i for i, v in enumerate(row) if v.strip() not in ('0', '0.0')]
if 0 not in keep:
    keep.insert(0, 0)
out = path.replace('_pathways.tsv', '_pathways_filtered.tsv')
with open(out, 'w') as f:
    f.write('\t'.join([hdrs[i] for i in keep]) + '\n')
    f.write('\t'.join([row[i] for i in keep]) + '\n')
print("Total: %d, Non-zero: %d -> %s" % (len(row) - 1, len(keep) - 1, out))
