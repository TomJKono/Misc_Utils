Misc_Utils
==========

A collection of Bash and Python scripts to perform various small tasks for DNA sequence data.

Contains the following scripts:
- **Count_Variants_Per_Contig.py**: Counts how many variants there are in each contig/chromosome in a VCF
- **Filter_VCF.py**: Apply arbitrary filters to a VCF file.
- **Genotype_Matrix_To_Fasta.py**: Convert a genotyping matrix to FASTA for input into [libsequence](http://molpopgen.github.io/libsequence/) tools. Because it assumes a fixed genotyping platform, it will remove monomorphic markers as well.
- **PLINK_to_NicholsonFST.py**: Convert from [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) file formats to those used as input for an *F*<sub>ST</sub> estimator developed by [Nicholson et al. 2002](http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00357/abstract) and implemented in the R package '[popgen](http://cran.r-project.org/web/packages/popgen/index.html)'
- **Percent_Similarity.py**: Calculate average percent similarity for SNPs listed in a VCF file.
- **Parallel_ms.py**: Splits a big ms simulation over multiple cores to make it run in less walltime.
- **Plot_SFS.R**: Plot site frequency spectra.
- **Remove_Monomorphic.py**: Remove monomorphic sites from a FASTA alignment.
- **Strip_BAM.sh**: Trim down a BAM file to just regions of interest. Requires [SAMTools](http://www.htslib.org).
- **VCF_MAF.py**: Counts the number of alternate and reference reads in a VCF. Useful only for BWC's BSA project (for now)
- **ms_FreqFilter.py**: Apply a 'discovery panel' to ms output. Used in [Fang et al. 2013](http://www.g3journal.org/content/3/11/1945.abstract) in G3 to simulate ascertainment for a genotyping platform.
- **transpose.sh**: Transpose a matrix.
