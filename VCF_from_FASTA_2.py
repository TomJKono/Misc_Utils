#!/usr/bin/env python
"""Updated version of VCF_from_FASTA.py that hopefully will simplify the
process of generating a VCF from Sanger reads aligned to a reference sequence.
Requires Biopython. Assumes the reference sequence is the first one in the
alignment. Takes one argument:
    1) FASTA multiple sequence alignment"""

import sys
import re
from datetime import date
try:
    from Bio import SeqIO
    fa_in = sys.argv[1]
except ImportError:
    sys.stderr.write('This script requires Biopython.\n')
    sys.exit(1)
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(2)


def extract_ref_coords(alignment, ref_idx=0):
    """Extract the name and the coordinates from the reference sequence of the
    alignment. This expects the sequence to be named in SAM/BED region
    format."""
    ref = alignment[ref_idx]
    chrom = ref.name.split(':')[0]
    start = ref.name.split(':')[1].split('-')[0]
    return (chrom, start)


def check_ref_gaps(alignment, ref_idx=0):
    """Check the reference sequence for end gaps. If we find any end gaps on
    the left, we will throw an error because we cannot accurately calculate the
    positions relative to reference. If we find them on the right, we will
    throw a warning, but continue anyway."""
    ref = alignment[ref_idx]
    refseq = str(ref.seq)
    left_gap = re.compile(r'^-+[ACGTMRWSYKVHDBNacgtmrwsykvhdbn]')
    right_gap = re.compile(r'[ACGTMRWSYKVHDBNacgtmrwsykvhdbn]-+$')
    if left_gap.match(refseq):
        sys.stderr.write(
            """Error:
Reference sequence has end-gaps on the left. This will cause calculated
positions to be incorrect. Please remove the end gaps and re-run this script.
""")
        sys.exit(10)
    if right_gap.search(refseq):
        sys.stderr.write(
            """Warning:
reference sequence has end-gaps on the right. This is not an error, but some of
the variants at the end of the alignment will not be placed on the reference
sequence. You may either remove the right end-gap or remove the variants from
the resulting VCF that occur in the end-gap.\n""")
    return


def extract_variants(alignment, ref_idx=0):
    """Extract the positions of SNPs and indels in the Sanger reads aligned to
    the reference sequence."""
    snp_pos = []
    indel_pos = []
    # First, convert the alignment to a list of lists, as opposed to a list
    # of SeqRecord objects
    raw_aln = [list(s.seq) for s in alignment]
    # Transpose it so that we iterate over columns of the alignment
    t_raw_aln = zip(*raw_aln)
    # Start iterating across columns and saving positions of variant sites. We
    # keep track of the reference base and only increment the position counter
    # for when we see a non-gap character in the reference sequence.
    offset = 0
    for aln_column in t_raw_aln:
        # First, get the states that exist at this position
        states = set(aln_column)
        # Discard any 'N' bases
        states.discard('N')
        # And get the ref state
        ref_state = aln_column[ref_idx]
        # Use the ref state to get the alternate states
        alt_states = states - set(ref_state)
        # If there is a gap in this position, then we will append it to the
        # list of indel positions
        if '-' in states:
            indel_pos.append((offset, ref_state, alt_states))
        # Then, discard the gap character to look for SNPs
        states.discard('-')
        # If the length of the states is greater than 1, then we have a SNP
        if len(states) > 1:
            # We will calculate the following:
            #   Number of non-missing alleles
            #   Minor allele count
            #   Minor allele frequency
            # The reference IS included in these calculations.
            non_missing = [
                base
                for base
                in aln_column
                if base != '-'
                or base != 'N']
            acs = [aln_column.count(x) for x in states]
            afs = [float(c)/len(non_missing) for c in acs]
            # Re-discard the gap character, just to be sure we do not count it
            # as an alternate state
            alt_states.discard('-')
            snp_pos.append(
                (offset, ref_state, alt_states, len(non_missing), min(acs),
                 min(afs)))
        # If the reference sequence is not a gap, then we increment our offset
        # counter.
        if ref_state != '-':
            offset += 1
    return (snp_pos, indel_pos)


def collapse_indels(indels):
    """Collapse indels by identifying runs of consecutive integers and merging
    those into a single entry."""
    # Sort the indel bases by their coordinate
    indel_srt = sorted(indels, key=lambda x: x[0])
    # Make a list to hold our aggregated indels
    agg_indel = []
    # We will now iterate over adjacent records - if they are consecutive, then
    # merge them. Else, break it and start a new record.
    curr_indel = []
    for ind, ind_adj in zip(indel_srt, indel_srt[1:]):
        # Unpack the alleles. It's a little silly, but we have to cast the set
        # to a list to subset it.
        curr_ref = ind[1]
        curr_alt = list(ind[2])[0]
        adj_ref = ind_adj[1]
        adj_alt = list(ind_adj[2])[0]
        if not curr_indel:
            curr_indel = [ind[0], curr_ref, curr_alt]
        # If the next position is not consecutive, append it and start over
        if ind_adj[0] - ind[0] > 1:
            agg_indel.append(curr_indel)
            curr_indel = [ind_adj[0], adj_ref, adj_alt]
        else:
            curr_indel[1] += adj_ref
            curr_indel[2] += adj_alt
    # The way we are iterating through the indel list means that we will always
    # leave off the last one. Append it after the loop finishes.
    agg_indel.append(curr_indel)
    return agg_indel


def adjust_indels(indels, alignment, ref_idx=0):
    """Adjust the indel positions so that they are offset by one, as required
    by the VCF spec. This is because the reported position must be the base
    *before* any insertion/deletion polymorphism is observed."""
    spec_indels = []
    # Remove the gaps from the reference sequnce for getting the reference base
    # of the indel
    ref_seq = ''.join([base for base in alignment[ref_idx].seq if base != '-'])
    for i in indels:
        spec_pos = i[0] - 1
        spec_ref = ref_seq[spec_pos]
        spec_indel = [spec_pos, spec_ref + i[1], spec_ref + i[2]]
        spec_indels.append(spec_indel)
    return spec_indels


def print_vcf(snp_var, ind_var, refseq, offset):
    """Print a VCF from the calculated positions of the variants."""
    # Define the VCF header
    today = date.today().strftime('%Y%m%d')
    vcf_header = """##fileformat=VCFv4.1
##fileDate={filedate}
##source=VCF_from_FASTA_2.py;Morrell Lab @ UMN
##INFO=<ID=MAC,Number=1,Type=Integer,Description="Minor allele count">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Variant is an INDEL">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO""".format(filedate=today)
    print(vcf_header)
    # Sort the SNPs and indels by their offset and print the records
    srt_variants = sorted(snp_var + ind_var, key=lambda x: x[0])
    for v in srt_variants:
        # Set the chromosome and position: add 1 to account for 1-based nature
        # of VCF
        v_chr = refseq
        v_pos = str(int(offset) + v[0] + 1)
        v_id = '.'
        v_qual = '40'
        v_filter = '.'
        # This is a bit of a hack, but if we have more than three fields, then
        # the variant type is a SNP
        if len(v) > 3:
            v_ref = v[1]
            # A bit ugly, but we have to cast the alt alleles from set to list
            v_alt = ','.join(list(v[2]))
            v_info = ';'.join([
                'MAC=' + str(v[4]),
                'MAF=' + str(v[5]),
                'NS=' + str(v[3]),
                'SNP'])
        else:
            # For indels, replace the gap characters with N
            v_ref = v[1].replace('-', 'N')
            v_alt = v[2].replace('-', 'N')
            v_info = 'INDEL'
        # Print the line
        print('\t'.join([
            v_chr,
            v_pos,
            v_id,
            v_ref,
            v_alt,
            v_qual,
            v_filter,
            v_info]))
    return


def main(fasta):
    """Main function."""
    # Store the alignment object as a list
    aln = list(SeqIO.parse(fasta, 'fasta'))
    # Extract the chromosome name and start position from the name of the
    # reference sequence.
    chrom, start = extract_ref_coords(aln)
    # We should check the reference sequence, too. If there are end-gaps on the
    # reference sequence, then we can't accurately calculate positions in the
    # alignment.
    check_ref_gaps(aln)
    snps, indels = extract_variants(aln)
    # Next, we want to collapse indels. We can find these by identifying runs
    # of consecutive integers in the list of indels. Some of the variants that
    # are in the list of indels are SNPs that occur within sequences that are
    # also part of a length polymorphism. We can just treat them as indels for
    # this routine.
    c_indels = collapse_indels(indels)
    # We also have to adjust the indels: the VCF spec requires that the
    # position of the indel be the base *before* the length polymorphism
    a_indels = adjust_indels(c_indels, aln)
    # Then, print the VCF!
    print_vcf(snps, a_indels, chrom, start)
    return


main(fa_in)
