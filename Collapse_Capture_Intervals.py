#!/usr/bin/env python
"""Collapse smaller consecutive intervals into larger ones. Allows a minimum
distance intervals are spaced for them to be consecutive."""

import sys

intervals = sys.argv[1]


def collapse_intervals(regions, dist, minsize=150):
    """Takes a list of tuples describing intervals, and collapses consecutive
    intervals into larger ones. If two intervals are 'dist' or fewer bp apart,
    they are considered to be in the same region. If an interval falls below a
    certain size, it is rejected."""
    new_reg = []
    #   First collapse consecutive intervals
    tmp_reg = ()
    for r in regions:
        if len(tmp_reg) == 0:
            tmp_reg += r
        else:
            if int(r[0]) - int(tmp_reg[-1]) > dist:
                new_reg.append((tmp_reg[0], tmp_reg[-1]))
                tmp_reg = ()
            else:
                tmp_reg += r

    #   Then filter by size
    new_reg_flt = [i for i in new_reg if int(i[1]) - int(i[0]) >= minsize]
    return(new_reg_flt)

chroms = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        if tmp[0] not in chroms:
            chroms[tmp[0]] = []
        chroms[tmp[0]].append((tmp[1], tmp[2]))

for c in sorted(chroms.keys()):
    captured = collapse_intervals(chroms[c], 50)
    for reg in captured:
        print '\t'.join([c, reg[0], reg[1]])
