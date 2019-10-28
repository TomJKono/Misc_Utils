#!/usr/bin/env python
"""Predicts the effect of a single nucleotide substitution on a protein
sequence. Takes a GFF, reference sequence, and VCF describing SNPs and returns
a file with predictions. Information returned includes the gene ID, transcript
ID, reside number of affected codon and whether or not the SNP is silent or
nonsynonymous. This will always report positions in 1-based coordinates, from
"start" to "end." This means that for reverse-strand features, the script will
count backwards.

Requires Biopython and gff_parse"""

#   Import some modules
import sys
import math
import gff_parse
#   Biopython modules. These handle sequence objects and translations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation


def read_sequence(refseq):
    """Reads the reference sequence and returns a dictionary of Seq obects."""
    handle = open(refseq, 'r')
    #   Start a parser that will iterate through the sequence file
    ref_parser = SeqIO.parse(handle, 'fasta')
    #   And then turn the parser into a sequence dictionary
    ref_dict = SeqIO.to_dict(ref_parser)
    return ref_dict


def read_gff(gff_file):
    """Reads and parses the GFF."""
    parsed_gff = gff_parse.GFFHandler()
    parsed_gff.gff_parse(gff_file)
    return parsed_gff


def build_cds_sequences(cds_chunks):
    """Takes a list of CDS annotations and sticks them together."""
    feature_parts = []
    for c in cds_chunks:
        #   We convert all the coordinates to 0-based, since they came out of
        #   the GFF as 1-based. SeqFeature.extract() will give an interval
        #   corresponding to [start, end), so we only have to subtract from
        #   the start.
        start = int(c.start) - 1
        end = int(c.end)
        if c.strand == '+':
            strand = 1
        elif c.strand == '-':
            strand = -1
        #   Then, we build a FeatureLocation object out of it, and give it the
        #   proper strandedness
        feature_sub = FeatureLocation(start, end)
        feature_parts.append(feature_sub)
    #   Check the strand, and sort accordingly. For forward strand features,
    #   we sort from low to high. For reverse strand features, we sort from
    #   high to low.
    if strand == 1:
        feature_parts.sort(key=lambda s: s.start)
    elif strand == -1:
        feature_parts.sort(key=lambda s: s.start, reverse=True)
    #   And then we concatenate them all together
    joint_feature = sum(feature_parts)
    joint_feature = SeqFeature(joint_feature, type='CDS', strand=strand)
    #   Return it!
    return joint_feature


def translate_codons(feature, sequence, alt, position):
    """Translates the two codon states and returns them."""
    #   First, extract the feature sequence from the main sequence
    #   Feature.extract() gives a SeqRecord, and we want just the Seq object
    #   inside it.
    cds_seq = feature.extract(sequence).seq
    seqlen = len(sequence)
    #   Then, we get the list of positions that are in the feature
    cds_positions = [x for x in range(0, seqlen) if x in feature]
    #   If the CDS is on the reverse strand, though, we need to reverse
    #   this list. We will return a flag, too, since we need to handle the
    #   ref and alt states accordingly, too
    rc = False
    if feature.strand == -1:
        cds_positions.reverse()
        #   We need to complement the alt state, too.
        alt = str(Seq(alt).complement())
        rc = True
    #   Sometimes, things slip through the filter and end up here. I'm not
    #   sure what causes this, but we need to trap this. In this case, we
    #   will return missing values for the effect predictions, since we can't
    #   accurately determine the effects.
    if position not in cds_positions:
        return(True, ['-', '-'], '-', '-')
    #   Then, which position in the CDS is the SNP?
    cds_base = cds_positions.index(position)
    #   Which codon is that?
    codon = int(math.floor(cds_base / 3))
    #   What are the positions covered by the codon?
    codon_bases = [
        cds_positions[codon*3],
        cds_positions[codon*3 + 1],
        cds_positions[codon*3 + 2]
        ]
    #   And what is the position of the SNP in the codon?
    snp_pos = codon_bases.index(position)
    #   Translate the reference CDS sequence
    ref_cds = cds_seq.translate()
    #   Then drop in the alternate base and re-translate
    alt_cds_seq = cds_seq
    #   make it a MutableSeq so that we can make substitutions
    alt_cds_seq = alt_cds_seq.tomutable()
    alt_cds_seq[cds_base] = alt
    alt_cds = alt_cds_seq.toseq().translate()
    states = [ref_cds[codon], alt_cds[codon]]
    #   Add 1 to make the offset 1-based
    residue_no = codon + 1
    return(rc, states, residue_no, snp_pos+1)


def usage():
    """Usage function."""
    message = """
Usage: SNP_Effect_Predictor.py [REF] [GFF] [VCF]

will predict the amino acid impact of SNPs listed in [VCF]. [REF] should be
a FASTA sequence, and [GFF] should reference the sequences in [REF].
Information returned is the SNP ID (Read from the VCF), the chromosomal
position, the gene ID (if genic), the transcript ID (if coding), whether the
alternate base represents a silent or nonsynonymous SNP, and the residue number
of the affected codon in 1-based coordinates.

Requires gff_parse.py from TomJKono's GitHub, and Biopython"""
    print message
    exit(1)


def main():
    """Main function."""
    if len(sys.argv[1:]) != 3:
        usage()
    #   Define arguments
    ref = sys.argv[1]
    gff = sys.argv[2]
    vcf = sys.argv[3]
    #   Parse the GFF information
    gff_data = read_gff(gff)
    #   And the reference sequence
    ref_dict = read_sequence(ref)
    #   Print out a header
    print '\t'.join(
        [
            'SNP_ID',
            'Chromosome',
            'Position',
            'Silent',
            'Transcript_ID',
            'Codon_Position',
            'Ref_Base',
            'Alt_Base',
            'AA1',
            'AA2',
            'CDS_Pos'
        ]
        )
    #   Start stepping through the VCF, and predicting the effects of each one
    with open(vcf, 'r') as f:
        for line in f:
            #   Skip comments and directives
            if line.startswith('#'):
                continue
            else:
                #   Separate the fields on tabs
                tmp = line.strip().split()
                #   We only want to save the first eight fields, as these are
                #   the ones that contain variant data (The others contain
                #   sample data.)
                chrom = tmp[0]
                #   Subtract 1 from the position to make it 0-based
                pos = int(tmp[1]) - 1
                snp_id = tmp[2]
                ref_base = tmp[3]
                alt_base = tmp[4]
                qual = tmp[5]
                flt = tmp[6]
                snp_info = tmp[7]
                #   We ask if either the ref or alt alleles involve more than
                #   one base. If this is the case, then we skip it, as we do
                #   not predict indel effects.
                if len(ref_base) > 1 or len(alt_base) > 1:
                    continue
                else:
                    #   Get all the CDS features that overlap the SNP
                    overlapping_cds = gff_data.overlapping_feature(
                        chrom,
                        pos,
                        feat_type='CDS')
                    #   Then, we get the sequence of the contig or chromosome
                    #   with the SNP
                    chrom_seq = ref_dict[chrom]
                    #   This is not the best approach, probably, but take the
                    #   first CDS feature that is overlapping
                    if overlapping_cds:
                        feat = overlapping_cds[0]
                        #   What transcript is it in?
                        transcript_id = gff_data.get_parents(chrom, feat.ID)
                        transcript_id = transcript_id[0].ID
                        #   Get the full CDS
                        full_cds = gff_data.get_siblings(
                            chrom,
                            feat.ID,
                            feat_type='CDS')
                        full_cds = build_cds_sequences(full_cds)
                        revcomp, states, aa_pos, codon_base = translate_codons(
                            full_cds,
                            chrom_seq,
                            alt_base,
                            pos)
                        #   Check if the SNP is synonymous or not
                        if len(set(states)) == 1:
                            silent = 'Yes'
                        else:
                            silent = 'No'
                        #   Check if we had to reverse complement or not
                        if revcomp:
                            ref_base = Seq(ref_base).complement()
                            alt_base = Seq(alt_base).complement()
                    else:
                        transcript_id = '-'
                        silent = 'Yes'
                        states = ['-', '-']
                        aa_pos = '-'
                        codon_base = '-'
                    #   Build the list of things to print
                    to_print = '\t'.join(
                        [
                            snp_id,
                            chrom,
                            str(pos + 1),
                            silent,
                            transcript_id,
                            str(codon_base),
                            str(ref_base),
                            str(alt_base),
                            states[0],
                            states[1],
                            str(aa_pos)
                        ]
                        )
                    print to_print
    return


main()
