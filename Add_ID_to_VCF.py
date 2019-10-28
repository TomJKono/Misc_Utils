#!/usr/bin/env python
"""Replace the ID column of a VCF with a sequential identifier, defined
internally. This is just to help keep track of variants throughout the project
that will get filtered for selected for various reasons. The IDs will be of the
form Prefix_Num, where Num is a zero-padded integer. Takes three arguments:
    1) VCF to process
    2) Prefix
    3) Number of zeroes to pad"""

import sys
import gzip

try:
    vcf = sys.argv[1]
    pref = sys.argv[2]
    pad = int(sys.argv[3])
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    exit(1)
except ValueError:
    sys.stderr.write('Please provide an integer for the pad value.\n')
    exit(2)


def main(v, pre, pval):
    """Main function."""
    varnum = 0
    with gzip.open(v, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                print(line.strip())
            else:
                tmp = line.strip().split('\t')
                snpid = pre + '_' + str(varnum).zfill(pval)
                newline = [tmp[0], tmp[1], snpid] + tmp[3:]
                print('\t'.join(newline))
                varnum += 1
    return


main(vcf, pref, pad)
