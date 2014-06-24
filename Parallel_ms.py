#!/usr/bin/env python

#   Takes three arguments:
#       1) The name of the MS script
#       2) The number of runs to do in parallel
#       3) The number of simulations per run
#       eg:
#           python Parallel_ms.py ms.sh 10 100
#       will run 'ms.sh' in 10 parallel jobs, each job doing 100 simulations

#   To spawn subprocesses
import subprocess
#   To take arguments
import sys
#   To deal with file paths
import os

#   This is the full path to the ms script
msscript = os.path.abspath(sys.argv[1])
#   The number of simulations to have per run
nsims = sys.argv[3]

#   A list for the handles to our subprocesses
processes = []
#   Go through and start our subprocesses
#   Keep track of them with the list
for i in range(1, int(sys.argv[2])+1):
    args = ['/bin/bash', msscript, str(i), nsims]
    p = subprocess.Popen(args)
    processes.append(p)

#   Now, we just wait for them to finish
for p in processes:
    p.wait()

#   Clean up a bit...
#   A handle to a new file that contains all the information we care about
handle = open('P_close', 'w')
#   Put the header in as expected
handle.write('pi\n')
#   For each of our temporary files, we read the lines, and put them at the end of the output
for i in range(1, int(sys.argv[2])+1):
    f = open('P_close_'+str(i), 'r')
    lines = f.readlines()
    handle.write(''.join(lines[1:]))
    f.close()
    #   And we delete the files
    os.remove('P_close_'+str(i))
    os.remove('prior_div_'+str(i))
    os.remove('prior_T_'+str(i))
    os.remove('prior2_'+str(i))

handle.close()
