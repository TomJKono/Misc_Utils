#!/usr/bin/env python
"""Randomly knock-down the quality of some proportion of FASTQ reads. Takes
one argument:
    1) FASTQ to process (gzipped)

Parameters for the quality alteration are given in the script as constants.
"""

import sys
import gzip
import numpy

# Define the FASTQ quality offset here
OFFSET = 33
# Define the proportion of reads to perturb as a probability
PROB_PERTURB = 0.2
# Define the number of bases to perturb. This will be counted in from the end
# of the read
NUM_BASES = 50
# Define the qualities to subtract and the probabilities that they should
# occur, if a read is selected to be perturbed
QUALS = [5, 10, 15, 20, 25]
QUAL_PROBS = [0.35, 0.25, 0.20, 0.1, 0.1]

# Iterate through the FASTQ, four lines at a time. We are making a possibly
# strong assumption here about the way the file is formatted. Each piece of
# information should be on one line.
with gzip.open(sys.argv[1], 'rt') as f:
    for seqname, nucleotides, comment, quals in zip(f, f, f, f):
        # First, decide if we are going to alter the quals
        perturb = numpy.random.binomial(1, PROB_PERTURB, size=1)[0]
        if perturb:
            sys.stderr.write(seqname.strip() + ' will be altered\n')
            # Decide how much we are going to perturb it
            amount = numpy.random.choice(QUALS, 1, p=QUAL_PROBS)[0]
            sys.stderr.write('Decreasing quality by ' + str(amount) + '\n')
            # Cast the quality string to decimals
            int_quals = [ord(c) - OFFSET for c in quals.strip()]
            # Perturb the qualities
            pert_quals = int_quals[:-NUM_BASES] + [q - amount for q in int_quals[-NUM_BASES:]]
            # Make sure that there aren't any negative qualities
            sane_quals = [q if q >= 0 else 0 for q in pert_quals]
            # Cast back to ascii
            new_quals = ''.join([chr(q+OFFSET) for q in sane_quals])
        else:
            new_quals = quals.strip()
        # Print out the new FASTQ
        toprint = [
            seqname.strip(),
            nucleotides.strip(),
            comment.strip(),
            new_quals.strip()]
        sys.stdout.write('\n'.join(toprint) + '\n')
