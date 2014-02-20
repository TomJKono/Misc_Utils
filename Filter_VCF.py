#!/usr/bin/env python

#   A script to apply various arbitrary filters to a VCF. The three filters
#   I have included here are genotype quality score, number of heterozygous
#   samples, number of samples with missing data, and read depth.
#   All threshholds can be adjusted by modifying the parameters at the top
#   of the script. This script was written with the VCF output from the
#   GATK version 2.8.1. The format may change in the future.
#   This script writes the filtered VCF lines to standard output

import sys
#   If variants have below a PHRED-scaled quality of 40,
#   we exclude them
quality_cutoff = 40
het_cutoff = 2
#   We exclude all missing data
missing_cutoff = 0
#   If we have low genotype confidence, then we also want to exclude the SNP
#   The genotype qualities are also stored as a PHRED-scaled probability
gt_cutoff = 40
n_gt_cutoff = 1
#   Our coverage cutoff is 5 reads per sample
per_sample_coverage_cutoff = 5
n_low_coverage_cutoff = 0
#   The number of samples
nsam = 1

#   Read the file in line-by-line
with open(sys.argv[1]) as f:
    for line in f:
        #   Skip the header lines - write them out without modification
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            tmp = line.strip().split('\t')
            #   we aren't confident in our ability to call ancestral state of
            #   indels
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue
            format = tmp[8]
            sample_information = tmp[9:]
            #   The genotype information is the first element of each sample
            #   info block in the VCF
            #   start counting up the number of hets, low quality genotypes and
            #   low coverage sites
            nhet = 0
            low_qual_gt = 0
            low_coverage = 0
            missing_data = 0
            for s in sample_information:
                #   For the GATK HaplotypeCaller, the per-sample information
                #   follows the form
                #   GT:AD:DP:GQ:PL
                info = s.split(':')
                if len(info) != 5:
                    break
                gt = info[0]
                #   We have to check for missing data first, because if it is
                #   missing, then the other fields are not filled in
                if '.' in gt:
                    missing_data += 1
                else:
                    dp = info[2]
                    gq = info[3]
                    #   Split on / (we have unphased genotypes) and count how
                    #   many hets we have
                    if len(set(gt.split('/'))) > 1:
                        nhet += 1
                    if dp == '.' or int(dp) < per_sample_coverage_cutoff:
                        low_coverage += 1
                    if gq == '.' or int(gq) < gt_cutoff:
                        low_qual_gt += 1
            #   The quality score is the sixth element
            if tmp[5] == '.' or float(tmp[5]) < quality_cutoff or nhet > het_cutoff or low_qual_gt > n_gt_cutoff or low_coverage > n_low_coverage_cutoff or missing_data > missing_cutoff:
                continue
            else:
                sys.stdout.write(line)
