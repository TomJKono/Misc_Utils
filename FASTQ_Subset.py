#!/usr/bin/env python
"""Subset a pair of FASTQ files based on a list of read names. Takes three
arguments:
    1) List of read names to keep
    2) R1 (gzipped)
    3) R2 (gzipped)
Will also take a predefined number of non-target reads, to keep the subset
looking like a "real" dataset. This is for a tutorial that uses a subset of
the genome.
"""

import sys
import gzip

# Set the number of non-target fragments here. The dataset is going to use a
# set of reads mapped to chromosome 19 of mm10. There are about 900,000
# fragments mapped to chr19 in the test dataset, so we will take 135,000
# non-target fragments (~15% of the library)
OFF_TARGET = 135000

# We are going to assume that the forward reads and the reverse reads differ
# by only the final character and that they are in the same order. This should
# be true for well-formed FASTQ files
TARGET_FRAGMENTS = set()
with open(sys.argv[1], 'r') as f:
    for line in f:
        fwd_read_name = '@' + line.strip() + '/1'
        TARGET_FRAGMENTS.add(fwd_read_name)

# Iterate through the R1 and R2 files to keep the target reads
r1 = gzip.open(sys.argv[2], 'rt')
r2 = gzip.open(sys.argv[3], 'rt')

r1_o = gzip.open(sys.argv[2].replace('.fastq.gz', '_Sub.fastq.gz'), 'wt')
r2_o = gzip.open(sys.argv[3].replace('.fastq.gz', '_Sub.fastq.gz'), 'wt')

KEPT_OFF_TARGET = 0
for r1_name, r1_nuc, r1_comment, r1_qual, r2_name, r2_nuc, r2_comment, r2_qual in zip(r1, r1, r1, r1, r2, r2, r2, r2):
    # first write a bunch of reads that don't map to chr19
    if KEPT_OFF_TARGET <= OFF_TARGET and r1_name.strip() not in TARGET_FRAGMENTS:
        curr_r1 = r1_name + r1_nuc + r1_comment + r1_qual
        curr_r2 = r2_name + r2_nuc + r2_comment + r2_qual
        r1_o.write(curr_r1)
        r2_o.write(curr_r2)
        KEPT_OFF_TARGET += 1
    else:
        # Then, check if R1 is in the set of target names
        if r1_name.strip() in TARGET_FRAGMENTS:
            curr_r1 = r1_name + r1_nuc + r1_comment + r1_qual
            curr_r2 = r2_name + r2_nuc + r2_comment + r2_qual
            r1_o.write(curr_r1)
            r2_o.write(curr_r2)
        else:
            continue

r1_o.flush()
r1_o.close()
r2_o.flush()
r2_o.close()
