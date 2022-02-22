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
print('Group'.ljust(14) + 'TB.Used'.ljust(8) + 'TB.Quota'.ljust(9) + '%TB'.ljust(6) + 'M.Files'.ljust(8) + 'F.Quota'.ljust(8) + '%Files'.ljust(8) + 'Scratch.TB'.ljust(11) + '%Scratch'.ljust(10))
for g in sorted(gp_list):
    gp_q_cmd = ['groupquota', '-f', '-g', g, '-c', '-U', 'T', '-H']
    # We need another command to get the scratch quotas, too.
    gp_q_cmd_scratch = ['groupquota', '-f', '-g', g, '-c', '-U', 'T', '-H', '-r']
    p = subprocess.Popen(gp_q_cmd, shell=False, stdout=subprocess.PIPE)
    quota = p.communicate()[0]
    quota = quota.decode('utf-8')
    # Run the scratch quota checker.
    p_scratch = subprocess.Popen(gp_q_cmd_scratch, shell=False, stdout=subprocess.PIPE)
    quota_scratch = p_scratch.communicate()[0]
    quota_scratch = quota_scratch.decode('utf-8')
    # Split it on commas and calculate the info
    group, tb_used, tb_lim, nfiles, flim = quota.split(',')
    group, scratch_used, scratch_lim, scratch_nfiles, scratch_flim = quota_scratch.split(',')
    nfiles = nfiles.split('.')[0]
    flim = flim.split('.')[0]
    # Replace the '000000' in N files with 'M'
    #   We have to do this dumb reverse-replace-reverse thing to make sure that
    #   it replaces the rightmost 000000 in the string. E.g., 10000000 -> 10M,
    #   rather than 10000000 -> 1M0.
    flim_prt = flim[::-1].replace('000000', 'M')[::-1]
    # Report files in millions, too
    nfiles_prt = str(round(float(nfiles) / 1000000, 2))
    try:
        perc_used = str(round((float(tb_used) / float(tb_lim)) * 100, 2))
    except ZeroDivisionError:
        perc_used = 'NA'
    try:
        perc_files = str(round((float(nfiles) / float(flim)) * 100, 2))
    except ZeroDivisionError:
        perc_files = 'NA'
    # And for scratch
    try:
        scratch_perc_used = str(round((float(scratch_used) / float(scratch_lim)) * 100, 2))
    except ZeroDivisionError:
        scratch_perc_used = 'NA'
    print(group.ljust(14) + tb_used.ljust(8) + tb_lim.ljust(9) + perc_used.ljust(6) + nfiles_prt.ljust(8) + flim_prt.ljust(8) + perc_files.ljust(8) + scratch_used.ljust(11) + scratch_perc_used.ljust(10))
