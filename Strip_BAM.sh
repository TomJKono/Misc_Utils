#!/bin/bash

#	A shell script to remove all but the specified contigs or chromosomes
#	from a given bam file. Requires samtools
#	arguments:
#		-r <regions file> a file containing the regions of interest, in a 
#			format that samtools understands:
#			contig:start-stop
#			contig2:start-stop
#			With one region per line
#		-o <output directory> the directory where the new BAM files will be 
#			written. This directory should be writable. This script creates
#			a sorted BAM file and an index.
#		-b <file> BAM file to operate on

#	We give a help message if no arguments are passed
SCRIPT=`basename $0`
if (($# == 0))
then
cat << EOF
Usage:
$SCRIPT -b <file> -r <regions file> -o <output dir>
Will remove all but the regions listed from the given BAM file, and convvert
it to BAM, sort it, then index it. All files will be written to the
<output dir> specified.
EOF
exit 1
fi

#	The path to samtools
#	If it's not in $PATH, then give it here
SAMTOOLS=samtools
#	A statement to check if samtools is installed, or correctly modified
command -v ${SAMTOOLS} > /dev/null 2>&1 || { echo >&2 "${SAMTOOLS} is not available! Is it properly set in the script?"; exit 1; }

#	 We will use getopts to parse the arguments
while getopts ":b:c:r:o:" OPT
do
	case $OPT in 
	r)
		#	Reading in the regions file to strip
		#	We only want to save the chromosomes (or contigs) since
		#	these will be used in the header
		#	We do this very ugly sed pipe sed for now, so that it leaves the word boundary
		#	delimiters off the very end of the string
		GREPREGIONS=`cut -d ':' -f 1 $OPTARG | tr '\n' ',' | sed -e 's/,/\\\b|/g' | sed -e 's/\\b|$//g'`
		#	This reads in the regions listed in the file, and changes
		#	all the newlines to spaces, which samtools will parse
		SAMREGIONS=`tr '\n' ' ' < $OPTARG`
		#	Formatting the regions for grep
		#GREPREGIONS=`echo $OPTARG | sed -e 's/,/\\\b|/g'`
		#	Formatting the regions for samtools
		#SAMREGIONS=`echo $OPTARG | sed -e 's/,/\ /g'`
		;;
	b)
		#	The full BAM file that we want to trim down
		BAMFILE=$OPTARG
		#	The filename that will output to
		#	replace the .bam extension with _trimmed.sam
		#	We use basename, since we want just the filename and not the
		#	directory that it is in
		OUTFILE=`basename ${OPTARG/.bam/_trimmed.sam}`
		;;
	o)
		OUTDIR=$OPTARG
		;;
	:)
		echo "Option -$OPTARG requires an argument!"
		exit 1
		;;
	\?)
		echo "-$OPTARG is not a valid argument!"
		exit 1
		;;
	esac
done

#	create this directory if it doesn't exist
mkdir -p ${OUTDIR}
#	These are lines we always want to save
SAVE='@PG\b|@RG\b|@HD\b|@CO\b|'
#	Output the header 
echo "Generating header and trimming ... "
$SAMTOOLS view -H $BAMFILE | grep -E "$SAVE$GREPREGIONS\\b" > ${OUTDIR}/${OUTFILE}
#	And trim the alignments
$SAMTOOLS view $BAMFILE ${SAMREGIONS} >> ${OUTDIR}/${OUTFILE}
echo "Done!"
#	cd into the output directory, so the paths won't get messed up
cd ${OUTDIR}
echo "cd into ${OUTDIR} ..."
#	Then, convert it to bam
echo "Converting ${OUTFILE} to ${OUTFILE/.sam/.bam} ... "
$SAMTOOLS view -bS ${OUTFILE} > ${OUTFILE/.sam/.bam}
echo "Done!"
#	Sort it
echo "Sorting ${OUTFILE/.sam/.bam} ... "
$SAMTOOLS sort ${OUTFILE/.sam/.bam} ${OUTFILE/_trimmed.sam/_sorted}
echo "Done!"
#	And index it
echo "Indexing ${OUTFILE/_trimmed.sam/_sorted.bam} ... "
$SAMTOOLS index ${OUTFILE/_trimmed.sam/_sorted.bam}
echo "Done!"

#	Clean up the temp files that we don't need anymore
#	Uncomment this if you don't want to see the 'trimmed' files in the output
#	directory.
######
######	WARNING
######	This command will remove *every* file in the output directory that
######	doesn't have 'sorted' in its name. DO NOT uncomment this if you want 
######	to output to a directory with files you want to save
######
#find ${OUTDIR} -type f -not -name '*sorted*' | xargs rm
