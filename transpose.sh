#!/bin/bash
#   an awk script to transpose a matrix
#   use at own risk, does not do dimension checking
#   Yanked from Awk and Sed book
exec awk '
NR == 1 {
    n = NF
    for (i = 1; i <= NF; i++)
        row[i] = $i
    next
}
{
    if (NF > n)
        n = NF
    for (i = 1; i <= NF; i++)
        row[i] = row[i] "\t" $i
}
END {
    for (i = 1; i <= n; i++)
        print row[i]
}' ${1+"$@"}
