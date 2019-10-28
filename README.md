Misc_Utils
==========

A collection of Bash and Python scripts to perform various small tasks for DNA sequence data.

Contains the following scripts:
- **Add_ID_to_VCF.py**: Add stable identifiers to the `ID` field of a VCF.
- **Count_Variants_Per_Contig.py**: Counts how many variants there are in each contig/chromosome in a VCF
- **Filter_VCF.py**: Apply arbitrary filters to a VCF file.
- **Genotype_Matrix_To_Fasta.py**: Convert a genotyping matrix to FASTA for input into [libsequence](http://molpopgen.github.io/libsequence/) tools. Because it assumes a fixed genotyping platform, it will remove monomorphic markers as well.
- **Mass_Job_Deletion.sh**: Delete all owned [MSI](https://www.msi.umn.edu/) jobs on the current server. Does not ask for confirmation, be careful.
- **PLINK_to_NicholsonFST.py**: Convert from [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) file formats to those used as input for an *F*<sub>ST</sub> estimator developed by [Nicholson et al. 2002](http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00357/abstract) and implemented in the R package '[popgen](http://cran.r-project.org/web/packages/popgen/index.html)'
- **Percent_Similarity.py**: Calculate average percent similarity for SNPs listed in a VCF file.
- **Parallel_ms.py**: Splits a big batch of [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) simulations over multiple cores to make it run in less walltime. 
- **Plot_SFS.R**: Plot site frequency spectra.
- **Remove_Monomorphic.py**: Remove monomorphic sites from a FASTA alignment.
- **SNP_Effect_Predictor.py**: Predicts silent/nonsynonymous SNPs in a VCF, given a GFF and a reference assembly. Requires gff_parse.py and [Biopython](http://biopython.org/). [SNPEff](http://snpeff.sourceforge.net/) does this, but SNP_Effect_Predictor.py was written to work with a genome with an incomplete assembly.
- **SRA_Fetch.sh**: Downloads .sra files from [NCBI's Short Read Archive](http://www.ncbi.nlm.nih.gov/sra) using [LFTP](http://lftp.yar.ru/). Can fetch based on Experiment number, Run number, Sample number, or Study number.
- **Strip_BAM.sh**: Trim down a BAM file to just regions of interest. Requires [SAMTools](http://www.htslib.org).
- **VCF_MAF.py**: Counts the number of alternate and reference reads in a VCF. Useful only for BWC's BSA project (for now)
- **VCF_To_Htable.py**: Translates a VCF into a Hudson-like polytable. Chokes on heterozygous sites.
- **gff_parse.py**: Python classes to try to make reading/fetching chunks of data from a GFF v3 file easier. Gets parent, child, and "sibling" features given a feature identifier.
- **ms_FreqFilter.py**: Apply a 'discovery panel' to [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) output. Used in [Fang et al. 2013](http://www.g3journal.org/content/3/11/1945.abstract) in G3 to simulate ascertainment for a genotyping platform.
- **transpose.sh**: Transpose a matrix.
