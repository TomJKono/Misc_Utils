#!/usr/bin/env python

#   A script to calculate the alt allele frequency in a VCF file
#   This is useful only for a BSA-Seq project
#   Uses the AD annotation from the GATK output

import sys

#   Start iterating through the file
with open(sys.argv[1], 'r') as f:
    for line in f:
        #   ignore header lines
        if line.startswith('##'):
            continue
        #   This defines how many samples in the VCF
        elif line.startswith('#CHROM'):
            samples = line.strip().split('\t')[9:]
            #   Create the list of sub-fields for each sample
            sample_sub_fields = ['_AltAlleleFreq', '_ReadDepth']
            #   And tack them together
            per_sample = []
            for s in samples:
                for sb in sample_sub_fields:
                    per_sample.append(s+sb)
            #   Append the header to the list of data to write
            print 'SNPPos\tRefAllele\tAltAllele\t' + '\t'.join(per_sample) + '\tNotes'
        else:
            tmp = line.strip().split('\t')
            #   Parse out the relevant information
            chromosome = tmp[0]
            bp_pos = tmp[1]
            ref = tmp[3]
            alt = tmp[4]
            filt = tmp[6]
            sample_info = tmp[9:]
            if 'LowQual' in filt:
                notes = 'Low Confidence Genotype'
            else:
                notes = ''
            format = tmp[8].split(':')
            sample_genotypes = [x.split(':') for x in sample_info]
            #   check if AD is not in the format field
            #   if not, then skip it
            if 'AD' not in format:
                notes = 'Missing Genotype Call'
                print '\t'.join([chromosome + ':' + bp_pos, ref, alt] + len(samples)*['NA', 'NA'] + [notes])
            else:
                notes = ''
                #   Which column is the AD?
                AD_pos = format.index('AD')
                #   For each sample...
                sample_variant_info = []
                for s in sample_info:
                    #   If it is totall missing
                    if s == './.':
                        AAF = 'NA'
                        total_depth = 'NA'
                    else:
                        gt = s.split(':')
                        #   If that sample has missing information
                        if gt[AD_pos] == '.':
                            AAF = 'NA'
                            total_depth = 'NA'
                        else:
                            #   Get the allele depths
                            AD = gt[AD_pos].split(',')
                            #   Convert to float
                            AD_flt = [float(i) for i in AD]
                            total_depth = sum(AD_flt)
                            if total_depth == 0:
                                AAF = '1.0'
                            else:
                                AAF = str(AD_flt[1]/total_depth)
                        sample_variant_info += [str(AAF), str(total_depth)]
                sample_list = [chromosome + ':' + bp_pos, ref, alt] + sample_variant_info + [notes]
                print '\t'.join(sample_list)
