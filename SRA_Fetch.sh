#!/bin/bash
#   Script to download all .sra files given an experiment number.
#   This is because the University of Minnesota Minnesota Supercomputing
#   Institute (MSI) does not allow remote access for the SRA Toolkit.
#   Requires the SRA Toolkit and lftp be installed.
#   Thomas Kono
#   Saint Paul, MN

set -e
set -u
set -o pipefail

USAGE="Usage:
$0 <Option> <Accession> [ -d DIR ] [ -v ] [ -h ]

Will fetch all .SRA files under the SRA accession number given by
<Accession>. If DIR is specified, then this script will create it if it does
not exist, and put everything in there. If DIR is not supplied, then the
current directory is used.

Pass -h to see available options."

HELP="Usage:
$0 <Option> <Accession> [ -d DIR ]

Available options:
Required:
    -e      Provided <Accession> is an Experiment number.
    -r      Provided <Accession> is a Run number.
    -p      Provided <Accession> is a Sample number.
    -s      Provided <Accession> is a Study number.

Optional:
    -d DIR  Output all .SRA files into DIR.

Switches:
    -v      Don't actually make directories and download, just print what is
            going to happen
    -h      Show this message and exit.
"

#   If there are no arguments passed to the script, drop the usage message and
#   exit.
if [ $# == 0 ]
    then echo "$USAGE"
    exit 1
fi

#   Parse the options
DRY_RUN="false"
DIR=$(pwd)
BASE_URL="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads"
while [[ $# > 0 ]]
    do
        FLAG="$1"
        case $FLAG in
            -e)
            ACC="$2"
            BASE_URL="${BASE_URL}/ByExp/sra"
            shift
            ;;
            -r)
            ACC="$2"
            BASE_URL="${BASE_URL}/ByRun/sra"
            shift
            ;;
            -p)
            ACC="$2"
            BASE_URL="${BASE_URL}/BySample/sra"
            shift
            ;;
            -s)
            ACC="$2"
            BASE_URL="${BASE_URL}/ByStudy/sra"
            shift
            ;;
            -d)
            DIR="$2"
            shift
            ;;
            -v)
            DRY_RUN="true"
            ;;
            -h)
            echo "$HELP"
            exit 2
            ;;
            *)
            echo "$USAGE"
            exit 1
            ;;
        esac
    shift
    done

#   Now that we have decided what type of accession number we have, we finish
#   building the rest of the URL
#   As of 2015-09-17, it follow this format:
#   ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/[Type]/sra/[SR_/ER_/DR_]
#       /[First 6 of Acc]/[Acc]/[Acc].sra
QUERY_URL="${BASE_URL}/${ACC:0:3}/${ACC:0:6}/${ACC}/"

#   if we are doing a dry-run, just print out these values and exit
if [[ "${DRY_RUN}" == "true" ]]
    then
        echo "The following operations will be performed:"
        echo "mkdir -p ${DIR}"
        echo "cd ${DIR}"
        echo "lftp -c mirror ${QUERY_URL} ."
        exit 3
    else
        #   Make the directory and cd into it
        mkdir -p ${DIR}
        cd ${DIR}
        #   get the contents with lftp
        lftp -c "mirror ${QUERY_URL} ."
fi
