#!/usr/bin/env python
"""Daft script to calculate average pairwise differences between all samples
listed in a VCF, accounting for the missing data among each comparison. The
diagonal is the number of called genotypes and the observed heterozygosity.
Takes two argument:
    1) VCF (gzipped)
    2) Output prefix"""

import gzip
import sys

OUT_PREFIX = sys.argv[2]

PI_OUT = OUT_PREFIX + '_Pi.csv'
NSITE_OUT = OUT_PREFIX + '_L.csv'


def make_square_mat(n, fill=0):
    """Return a square matrix of dimension n, populated with the specified
    fill value."""
    m = [
            [fill for i in range(0, n)]
            for j
            in range(0, n)
        ]
    return m


def pairwise_diff(g1, g2, auto=False):
    """Given two diploid genotype calls as strings, return the average pairwise
    differences. If auto=True, then we treat the genotypes as coming from the
    same individual. It will effectively return 1 for heterozygous sites and 0
    for homozygous sites."""
    # If either of the genotype calls are missing, return a value that will
    # let us skip the comparison without affecting the denominator.
    if g1 == './.' or g2 == './.':
        return None
    # Next easiest case, check if auto=True. If the alleles are identical,
    # return 0, else return 1
    if auto:
        g1_alleles = g1.split('/')
        if g1_alleles[0] == g1_alleles[1]:
            return 0.0
        else:
            return 1.0
    # We are assuming that the calls are all 0/1.
    ref = g1.count('0') + g2.count('0')
    alt = g1.count('1') + g2.count('1')
    # This is a bit lazy - we will hard-code in pi values for the allele
    # combinations we can have
    if ref == 4 or alt == 4:
        return 0.0
    elif ref == 1 or alt == 1:
        return 0.5
    elif ref == 2:
        return 2/3
    else:
        return None


# start iterating through the VCF
nsites_proc = 0
with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            hdr = line.strip().split('\t')
            samples = hdr[9:]
            nsamp = len(samples)
            # Make square matrices for both the number of sites and the number
            # of pairwise differences.
            l_mat = make_square_mat(nsamp)
            d_mat = make_square_mat(nsamp)
        else:
            if nsites_proc % 1000 == 0:
                sys.stderr.write('Processed ' + str(nsites_proc) + ' sites.\n')
            tmp = line.strip().split('\t')
            # Pass a replacement statement over the genotypes to replace them
            # with values that can be easily processed
            gts = [
                './.'
                if x.split(':')[0] == '.'
                else x.split(':')[0]
                for x
                in tmp[9:]
                ]
            # Iterate through all pairs of the genotypes and calculate the
            # average pairwise difference
            for i, geno1 in enumerate(gts):
                for j, geno2 in enumerate(gts):
                    # Calculate the pairwise difference
                    if i == j:
                        apd = pairwise_diff(geno1, geno2, auto=True)
                    else:
                        apd = pairwise_diff(geno1, geno2, auto=False)
                    # If the pairwise difference value is defined, then we have
                    # a case where we can increment the nsites and diff mat.
                    # Remember that 0.0 evaluates to boolean False. Ugh.
                    if apd is not None:
                        l_mat[i][j] += 1
                        d_mat[i][j] += apd
            nsites_proc += 1

# Now we want to print out the matrix
handle = open(NSITE_OUT, 'w')
handle.write(','.join([''] + samples) + '\n')
for s, row in zip(samples, l_mat):
    towrite = ','.join([s] + [str(x) for x in row]) + '\n'
    handle.write(towrite)
handle.close()

handle = open(PI_OUT, 'w')
handle.write(','.join([''] + samples) + '\n')
for s, row in zip(samples, d_mat):
    towrite = ','.join([s] + [str(x) for x in row]) + '\n'
    handle.write(towrite)
handle.close()
