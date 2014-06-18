#!/usr/bin/env python

#   A script to calculate the minor allele frequency in a VCF file
#   This is useful only for a single-sample VCF, for a BSA-Seq project
#   Uses the AD annotation from the GATK output

import sys

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            #   Parse out the relevant information
            chromosome = tmp[0]
            bp_pos = tmp[1]
            ref = tmp[3]
            alt = tmp[4]
            #   Insertion? Deletion? SNP?
            if len(ref) < len(alt):
                notes = 'Insertion WRT ref'
            elif len(ref) > len(alt):
                notes = 'Deletion WRT ref'
            else:
                notes = 'SNP'
            format = tmp[8].split(':')
            sample_info = tmp[9].split(':')
            #   check if AD is not in the format field
            #   if not, then skip it
            if 'AD' not in format:
                continue
            else:
                #   Which column is the AD?
                AD_pos = format.index('AD')
                #   Get the allele depths
                AD = sample_info[AD_pos].split(',')
                #   Convert to float
                AD_flt = [float(i) for i in AD]
                total_depth = sum(AD_flt)
                #   set a cutoff. Less than 10 reads -> no call
                if total_depth < 10:
                    continue
                else:
                    MAF = str(min(AD_flt)/total_depth)
                    #   Build the string to print out
                    to_print = '\t'.join([chromosome + ':' + bp_pos, ref, alt, MAF, str(int(total_depth)), notes])
                    print to_print
