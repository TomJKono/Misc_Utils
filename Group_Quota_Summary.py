#!/usr/bin/env python
"""Silly script to get a summary of the group quota for all of the groups that
a user is part of. Excludes software license groups."""

import subprocess

# Skip these groups
IGNORE = ['rnammer', 'signalp4', 'genemark']

# Command to list all groups
group_cmd = ['id', '-Gn']

# Run the group command and get the list
p = subprocess.Popen(group_cmd, shell=False, stdout=subprocess.PIPE)
gp_list = p.communicate()[0]
gp_list = gp_list.decode('utf-8')
gp_list = [g for g in gp_list.split() if g not in IGNORE]

# For each group, run the groupquota command with the -f flag to use cached
# data
print('Group'.ljust(16) + 'TB.Used'.ljust(10) + 'TB.Quota'.ljust(10) + '% Used'.ljust(10) + 'N.Files'.ljust(10) + 'File.Quota'.ljust(12) + '% Files'.ljust(10))
for g in sorted(gp_list):
    gp_q_cmd = ['groupquota', '-f', '-g', g, '-c', '-U', 'T', '-H']
    p = subprocess.Popen(gp_q_cmd, shell=False, stdout=subprocess.PIPE)
    quota = p.communicate()[0]
    quota = quota.decode('utf-8')
    # Split it on commas and calculate the info
    group, used, galaxy, tot, lim, nfiles, flim = quota.split(',')
    nfiles = nfiles.split('.')[0]
    flim = flim.split('.')[0]
    try:
        perc_used = round((float(tot) / float(lim)) * 100, 2)
    except ZeroDivisionError:
        perc_used = 'NA'
    try:
        perc_files = round((float(nfiles) / float(flim)) * 100, 2)
    except ZeroDivisionError:
        perc_files = 'NA'
    p_used = round(float(used), 2)
    print(group.ljust(16) + str(p_used).ljust(10) + lim.ljust(10) + str(perc_used).ljust(10) + nfiles.ljust(10) + flim.ljust(12) + str(perc_files).ljust(10))
