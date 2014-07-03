Misc_Utils
==========

A collection of Bash and Python scripts to perform various small tasks for DNA sequence data.

Contains the following scripts:
- Count_Variants_Per_Contig.py: Counts how many variants there are in each contig/chromosome in a VCF
- Filter_VCF.py: A python script to apply arbitrary filters to a VCF file.
- Percent_Similarity.py: A python script that calculates average percent similarity for SNPs listed in a VCF file.
- Parallel_ms.py: Splits a big ms simulation over multiple cores to make it run in less walltime.
- Plot_SFS.R: An R script to plot site frequency spectra.
- Remove_Monomorphic.py: A Python script to remove monomorphic sites from a FASTA alignment.
- Strip_BAM.sh: A bash script that trims down a BAM file to just regions of interest. Requires SAMTools.
- VCF_MAF.py: Counts the number of alternate and reference reads in a VCF. Useful only for BWC's BSA project (for now)
- ms_FreqFilter.py: Apply a 'discovery panel' to ms output. Used in Fang et al. 2013 in G3 to simulate ascertainment for a genotyping platform.
