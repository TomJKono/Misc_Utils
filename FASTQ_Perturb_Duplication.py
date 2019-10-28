#!/usr/bin/env python
"""Randomly add sequence duplication to a FASTQ.
    1) R1 FASTQ to process (gzipped)
    2) R2 FASTQ to process (gzipped)

Parameters for the duplication are defined as constants within the script.
"""

import sys
import gzip
import numpy

# Define the probability that a sequence will be duplicated
PROB_DUP = 0.0001
# Define the sequence duplication levels that are to be expected
DUP_LEVELS = [10, 50, 100, 500, 10000, 15000]
# And their probabilities
DUP_PROBS = [0.1, 0.1, 0.3, 0.3, 0.1, 0.1]

# Iterate through the two files
r1 = gzip.open(sys.argv[1], 'rt')
r2 = gzip.open(sys.argv[2], 'rt')

# Make output handles based on the filenames of the inputs
r1_o = gzip.open(sys.argv[1].replace('.fastq.gz', '_Dup.fastq.gz'), 'wt')
r2_o = gzip.open(sys.argv[2].replace('.fastq.gz', '_Dup.fastq.gz'), 'wt')
# Then, iterate through the reads in both
for r1_name, r1_nuc, r1_comment, r1_qual, r2_name, r2_nuc, r2_comment, r2_qual in zip(r1, r1, r1, r1, r2, r2, r2, r2):
    # decide if we are going to duplicate
    dup = numpy.random.binomial(1, PROB_DUP, size=1)[0]
    if dup:
        sys.stderr.write(r1_name.strip() + ' and ' + r2_name.strip() + ' will be duplicated\n')
        amount = numpy.random.choice(DUP_LEVELS, 1, p=DUP_PROBS)[0]
        sys.stderr.write('Duplication level: ' + str(amount) + '\n')
        numdup = 0
        while numdup <= amount:
            try:
                curr_r1 = r1_name + r1_nuc + r1_comment + r1_qual
                curr_r2 = r2_name + r2_nuc + r2_comment + r2_qual
                r1_o.write(curr_r1)
                r2_o.write(curr_r2)
                new_r1_name = next(r1)
                new_r1_nuc = next(r1)
                new_r1_com = next(r1)
                new_r1_qual = next(r1)
                new_r2_name = next(r2)
                new_r2_nuc = next(r2)
                new_r2_com = next(r2)
                new_r2_qual = next(r2)
                new_r1 = new_r1_name + r1_nuc + r1_comment + r1_qual
                new_r2 = new_r2_name + r2_nuc + r2_comment + r2_qual
                r1_o.write(new_r1)
                r2_o.write(new_r2)
                numdup += 1
            except StopIteration:
                numdup = 999999999
    else:
        new_r1 = r1_name + r1_nuc + r1_comment + r1_qual
        new_r2 = r2_name + r2_nuc + r2_comment + r2_qual
        r1_o.write(new_r1)
        r2_o.write(new_r2)

r1_o.flush()
r2_o.flush()
r1_o.close()
r2_o.close()
