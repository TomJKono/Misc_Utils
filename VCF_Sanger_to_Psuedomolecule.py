#!/usr/bin/env python

#   Creates a really crude VCF file given a FASTA alignment
#   Assumes that the reference sequence is the first sequence in the alignment

#   To take arguments
import sys
#   To handle ranges and getting elements
import itertools
import operator
import re
#   To read FASTA
from Bio import SeqIO

alignment = sys.argv[1]
sequences = list(SeqIO.parse(alignment, 'fasta'))
#   The reference is the first sequence in the alignment
reference = sequences[0]
#   And the name of the 'chromosome' is the contig name
string_name = reference.name
chrom_name = re.split(':|-', string_name)[0]
start_pos = int(re.split(':|-', string_name)[1]) + 2 # Add 2 because to convert from 0-based BED to 1-based VCF
#sys.exit(str(start_pos))

#   The VCF header
vcf_header="""##fileformat=VCFv4.1
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Variant is an INDEL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""

#   We then strip out just the sequence, leaving names etc. behind
#   Convert to a list, so we can save column information, too
raw_seqs = [list(x.seq) for x in sequences]
#   Transpose the matrix, to iterate over columns instead of rows
columns = zip(*raw_seqs)

#   Now, we actually find indels and SNP positions
indel_positions = []
snp_positions = []
#   For each column
for index, base in enumerate(columns):
    #   Convert to a set, which gives us how many states in each column
    tmp = set(base)
    #   Then we throw out Ns
    tmp.discard('N')
    #   this is a little ugly, but it works
    if '-' in tmp:
        indel_positions.append(index)
    #   Now, we discard the gaps, and count up SNPs
    tmp.discard('-')
    #   If there is more than one base, after throwing out gaps and Ns,
    #   then we have a SNP
    if len(tmp) > 1:
        snp_positions.append(index)

#   Now, we have a list of positions that are covered by SNPs and indels
#   We go through the indels and group them by consecutive ranges
#   This identifies runs of consecutive integers, and puts them into separate lists
#   for group_key, group in (0, pos1) (1, pos2) (2, pos3) ...
#   group by the value of (pos - index)
#   if they are the same, then the numbers are consecutive
indels = []
for key, group in itertools.groupby(enumerate(indel_positions), lambda (i,x):i-x):
    #   save the second value, which is the actual group
    indels.append(map(operator.itemgetter(1), group))

#   Now that we have a list of actual indels, we iterate over that, and save
#   the states that aren't gaps. We also save the reference state
indelfilename = alignment.replace('.fasta', '_INDELs.vcf')
handle = open(indelfilename, 'w')
handle.write(vcf_header)
for i in indels:
    #   First, we have to get the base before the indel
    #   to save it in the VCF
    #   Since this isn't in the indel positions, it should never be a gap
    before_indel = i[0]-1
    before_indel_base = reference[before_indel]
    #   Then, we get the reference and alternate states
    ref_allele = []
    alt_allele = []
    for position in i:
        ref_allele.append(reference[position])
        #   Get the alignment column, throw out Ns
        states = set(columns[position])
        states.discard('N')
        #   Then get the one that isn't the reference base
        if reference[position] == '-':
            alt_base = [s for s in states if s != reference[position]][0]
            alt_allele.append(alt_base)
        else:
            alt_allele.append('-')
    #   Now we have the reference allels and the alternate alleles
    #   We take out the gaps and find out which is 'longer'
    ref_allele = [x for x in ref_allele if x != '-']
    alt_allele = [x for x in alt_allele if x != '-']
    #   If the reference is longer, then we fill it with Ns, and truncate alt
    if len(ref_allele) > len(alt_allele):
        reference_state = before_indel_base + 'N'*len(i)
        alternate_state = before_indel_base
    else:
        #   Otherwise, truncate ref, and fill the alt with N
        reference_state = before_indel_base
        alternate_state = before_indel_base + 'N'*len(i)
    entry = '\t'.join([chrom_name, str(before_indel + start_pos), '.', reference_state, alternate_state, '.', '.', 'INDEL'])
    handle.write(entry + '\n')
handle.close()

#   Now, we make the SNP VCF files
snpfilename = alignment.replace('.fasta', '_SNPs.vcf')
handle = open(snpfilename, 'w')
handle.write(vcf_header)
for position in snp_positions:
    states = set(columns[position])
    states.discard('N')
    states.discard('-')
    reference_allele = reference[position]
    alternate_allele = [x for x in states if x != reference_allele][0]
    #   if it is involved in an indel, then we don't want it
    if reference_allele == '-' or alternate_allele == '-':
        continue
    entry = '\t'.join([chrom_name, str(position + start_pos), '.', reference_allele, alternate_allele, '.', '.', '.'])
    handle.write(entry + '\n')
handle.close()
