#!/usr/bin/env python3

import sys
import numpy as np
from collections import defaultdict

ld_file = sys.argv[1]
bin_size = 250
max_dist = 500000

bin_r2 = defaultdict(list)

with open(ld_file, "r") as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split()
        if len(parts) < 7 or parts[0] != parts[3]:
            continue
        try:
            bp1 = int(parts[1])
            bp2 = int(parts[4])
            r2 = float(parts[6])
        except ValueError:
            continue

        dist = abs(bp2 - bp1)
        if dist > max_dist:
            continue
        bin_index = dist // bin_size
        bin_r2[bin_index].append(r2)

# Output summary to TSV
with open("ld_decay_summary.tsv", "w") as out:
    out.write("Distance\tMean_r2\tMedian_r2\tCount\n")
    for bin_index in sorted(bin_r2):
        values = bin_r2[bin_index]
        mid = (bin_index + 0.5) * bin_size
        out.write(f"{mid}\t{np.mean(values):.4f}\t{np.median(values):.4f}\t{len(values)}\n")
