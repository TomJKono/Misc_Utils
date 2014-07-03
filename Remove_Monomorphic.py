#!/usr/bin/env python

#   A script to remove monomorphic sites from a FASTA alignment.
#   This is probably of limited use, but it still comes up every once
#   in a while.
#   Takes one argument:
#       1) The FASTA alignment
#   Takes one optional flag:
#       -i omit gaps (indels)
#   Outputs un-wrapped to stdout

import sys
import argparse

parser = argparse.ArgumentParser()
#   Add the alignment as a required argument
parser.add_argument('alignment', help='The FASTA alignment file')
#   And add the optional indel flag
parser.add_argument('-i', '--ignore-indels', help='Ignore gaps or indels', action='store_true')
args = parser.parse_args()


#   empty lists to hold the sequences and names from the FASTA file
names = []
sequences = []
with open(args.alignment, 'r') as f:
    for index, line in enumerate(f):
        #   if the line starts with >, then it is a new sequence
        if line.startswith('>'):
            #   If we are at the first line of the file, then we don't have
            #   any previously-parsed sequence to add
            if index != 0:
                sequences.append(seq.upper())
            seq = ''
            name = line.strip()
            names.append(name)
        else:
            #   Start accumulating the sequence
            #   We have to do it this way because the sequence may be wrapped
            #   which would throw off the line-by-line parsing
            seq += line.strip()
    #   Append the final seq of sequences
    sequences.append(seq.upper())

#   Run a quick check to see if the sequences are all the same length
first_seq = len(sequences[0])
#   True if all sequences are the same length, False otherwise
if not all(len(i) == first_seq for i in sequences):
    sys.stderr.write('Error! Your sequences are not all the same length! Are they aligned?\n')
    exit(1)

#   Next, we transpose the list of sequences to start removing columns
sequences_trans = zip(*sequences)
polymorphic = []
for column in sequences_trans:
    #   we use a set() type since it automatically removes duplicate items
    alleles = set(column)
    #   If we only have one base at that site, then we omit it
    if len(alleles) == 1:
        continue
    else:
        #   If there's a gap character '-' in this set, we skip over it
        #   but only if the -i flag is given
        if args.ignore_indels:
            if '-' in alleles:
                continue
            else:
                polymorphic.append(column)
        else:
            polymorphic.append(column)

#   Transpose it again to put them back in indivudal-rows
polymorphic_trans = zip(*polymorphic)
#   then print it out
for seqname, seq in zip(names, polymorphic_trans):
    print seqname
    print ''.join(list(seq))
