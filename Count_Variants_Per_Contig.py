#!/usr/bin/env python

#   A script to count the number of variants per contig for the Morex genome
#   assembly. This doesn't really make sense to do for a complete reference

import sys

#   create a dictionary
#   this won't matter because we don't have order anyway
contigs = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            chrom = tmp[0]
            if chrom not in contigs:
                #   If it's a contig we are seeing for the first time, then we
                #   add a new value to the dictionary, starting at 1 since this
                #   is a real count
                contigs[chrom] = 1
            else:
                #   Otherwise, we have seen it already, and we increment it
                contigs[chrom] += 1

#   Print it out
#   use iteritems() since it creates an iterator rather than a list
for ctg, count in contigs.iteritems():
    print ctg + '\t' + str(count)
