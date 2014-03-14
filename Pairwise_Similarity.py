#!/usr/bin/env python

#   A script to identify how many SNPs are private and shared, on average
#   between all pairwise comparisons in the samples in a VCF file

import sys

#   A function to calculate pairwise similarity
def pairwise(l1, l2):
    if len(l1) != len(l2):
        return(None)
    comparisons = 0
    similarity = 0
    for i in zip(l1, l2):
        if 'N' in i:
            continue
        else:
            comparisons += 1
            if i[0] == i[1]:
                similarity += 1
            else:
                continue
    return(similarity/float(comparisons))


#   A function to calculate the average
def average(x):
    return(sum(x)/float(len(x)))


vcf_data = []
#   Read in the VCF and store the data
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            tmp = line.strip().split('\t')
            #   This is the header line, and the sample information is
            #   on the 10th element until the end
            samples = tmp[9:]
            #   We will append each sample name to the data list
            #   but we will do it as a list with one element, so we can
            #   grow it as we read through the data
            for s in samples:
                vcf_data.append([s])
        else:
            #   We start reading through the variant call lines
            tmp = line.strip().split('\t')
            #   isolate the genotype data
            sample_information = tmp[9:]
            #   Append the genotype calls, 1 or 0 to the data list
            for index, s in enumerate(sample_information):
                genotype = s.split(':')[0]
                calls = genotype.split('/')
                #   Get the alleles
                alleles = set(calls)
                #   If we have a heterozygote, we append a missing data value
                if len(alleles) > 1 or '.' in alleles:
                    vcf_data[index].append('N')
                else:
                    #   get the proper call, 0 or 1
                    base = list(alleles)[0]
                    #   Append it!
                    vcf_data[index].append(base)

pairwise_similarity = []
#   Now, we go through and make all pairwise comparisons, counting up the
#   number of common and different alleles
for index, i in enumerate(vcf_data[:-1]):
    for j in vcf_data[index+1:]:
        #   We slice off the first element, since that is the sample name
        sim = pairwise(i[1:], j[1:])
        pairwise_similarity.append(sim)
        print i[0] + '-' + j[0] + ': ' + str(sim)

print 'Average: ' + str(average(pairwise_similarity))
