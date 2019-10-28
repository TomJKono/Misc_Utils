#!/usr/bin/env python
"""Calculates average nuceotide dieversity (Tajima's Pi) over a region using a
VCF and BED file. Assumes highly inbred samples, so genotypes are treated as
haploid calls, and heteroygous calls are treated as missing data. Note this
designed to be run over a single contig/chromosome, as it uses integers as
keys for identifying the number of unique variants."""

import sys


def n_choose_r(n, r):
    """nCr implemented recursively. Not the *MOST* efficient, but better than
    computing factorials."""
    if r == 1:
        return n
    elif r == n:
        return 1
    elif r == 0 or n == 0:
        return 0
    else:
        return n_choose_r(n-1, r-1) + n_choose_r(n-1, r)


def pairwise_diversity(calls):
    """Calculate pairwise diversity given a list of genotype calls. Returns 0
    for a monomorphic site, or sites with only one non-missing call."""
    #   Count up the number of reference and alternate genotypes.
    if 0 in calls:
        ref_count = calls.count(0)
    else:
        return 0
    if 1 in calls:
        alt_count = calls.count(1)
    else:
        return 0
    #   This sample size will change depending on how many non-missing genotypes
    #   there are.
    total_count = ref_count + alt_count
    #   Calculate up the similarities based on the number of reference and
    #   alternative genotypes. Calculate the number of pairwise comparisons
    #   that were made.
    ref_sim = n_choose_r(ref_count, 2)
    alt_sim = n_choose_r(alt_count, 2)
    total_comp = n_choose_r(total_count, 2)
    #   Then pairwise diversity is 1-[(ref_sim + alt_sim)/total_comp]
    return 1 - ((ref_sim + alt_sim) / float(total_comp))


def read_vcf(vcf):
    """Read a VCF, and calculate pairwise diversity for each site. Returns a
    a dictionary with genomic coordinate as key and pairwise diversity as the
    value."""
    vcfdata = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                pos = int(tmp[1])
                sample_info = tmp[9:]
                #   Get the genotype calls from the sample info fields
                gt = [t.split(':')[0] for t in sample_info]
                hap_calls = []
                for g in gt:
                    if g == '0/0':
                        hap_calls.append(0)
                    elif g == '1/1':
                        hap_calls.append(1)
                    else:
                        hap_calls.append('NA')
                site_pi = pairwise_diversity(hap_calls)
                vcfdata[pos] = site_pi
    return vcfdata


def calc_div_bed(vcfdata, bed):
    """Read a BED file, and for each interval, calculate the average pairwise
    diversity in it. Prints the interval start, end, and the average pairwise
    diversity."""
    avg_pairwise_divs = []
    with open(bed, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            #   Add 1 to the end, for 0-based BED coords
            start = int(tmp[1])
            end = int(tmp[2])
            #   Some intervals are of 0 length - return 0 for these.
            if start == end:
                continue
            #   Get the pi values that are between the current BED interval
            in_interval = [
                vcfdata[k]
                for k in vcfdata.keys()
                if (k >= start and k <= end)]
            #   Next, find the average pairwise diversity. We will just sum the
            #   diversities of the individual variants, then divide by the
            #   number of sites in the interval.
            pair_div = sum(in_interval) / float(end - start - len(in_interval))
            #   Then, append it to the list
            avg_pairwise_divs.append((str(start), str(end), str(pair_div)))
    return avg_pairwise_divs


def main(vcf, bed):
    """Main function."""
    v = read_vcf(vcf)
    d = calc_div_bed(v, bed)
    for i in d:
        print '\t'.join(i)
    return


main(sys.argv[1], sys.argv[2])
