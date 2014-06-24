#!/usr/bin/env python

#####
#   Importing modules
#####
#   To grab stdin
import sys
#   To get arguments
import argparse

#####
#   Functions
#####
def MAF(SNP):
    """This takes SNP genotypes in a list and returns the MAF."""
    return(min(SNP.count('1'), SNP.count('0')))

def mark_columns(panel, maf_floor):
    """This function takes the panel and the MAF floor, and marks columns
    for deletion out of the simulation. It returns a list of Y/N for
    columns that should be omitted."""
    #   We want to transpose the matrix, so each row is a SNP
    transposed_panel = zip(*panel)
    #   Iterate through the transposed matrix, build a list for positions
    #   that will be discarded because of MAF filter
    marks = []
    for each in transposed_panel:
        #   If the MAF is less than the floor
        if MAF(each) < maf_floor:
            #   Mark the position with a 'Y' -> this column will be trimmed
            marks.append('Y')
        else:
            #   If not, then mark with 'N' -> this column will be kept
            marks.append('N')
    return(marks)

def trim_simulations(simulation, panel_size, maf_floor):
    """This function actually trims columns out of the simulation. It also 
    calculates the new number of segregating sites, and prints out the updated
    list of positions."""
    #   Empty list for genotype matrix
    genotype_matrix = []
    #   If we have a block containing simulation data, then we iterate
    #   through it, and process accordingly
    for line in simulation:
        #   Save the positions list
        if line.startswith('pos'):
            #   We omit the first field, which would be 'positions:', 
            #   since we are only interested in the numbers
            positions = line.strip().split()[1:]
        #   If the line starts with 0 or 1, then we are in the actual
        #   SNP data from the simulation, we put that aside for processing
        if line.startswith('0') or line.startswith('1'):
            genotype_matrix.append(line.strip())
    #   Now, we get the columns that are to be filtered out
    marks = mark_columns(genotype_matrix[0:panel_size], maf_floor)
    #   The number of segregating sites is the number of times 'N' appears
    segsites = 'segsites: ' + str(marks.count('N'))
    #   We have to transpose the genotypes to iterate over columns
    transposed_genotypes = zip(*genotype_matrix)
    #   Then we traverse the whole thing, and only save those with an 'N'
    #   We start off the positions list with the correct label
    filtered_positions = ['positions:']
    filtered_genotypes = []
    for filt, pos, gen in zip(marks, positions, transposed_genotypes):
        if filt == 'N':
            filtered_positions.append(pos)
            filtered_genotypes.append(gen)
        else:
            continue
    #   Then, we want to put the genotype list back in the correct
    #   orientation
    filtered_genotypes = zip(*filtered_genotypes)
    #   And put it back into strings
    filtered_genotypes = [''.join(each) for each in filtered_genotypes]
    #   Start building the new simulation
    #   it will look like
    #
    #   segsites: nsegsites
    #   positions: 0.1 0.2 0.3 0.4 ...
    #   001101...
    #   100100...
    #   110010...
    #   ...
    #   
    #   //
    new_sim = ['//', segsites, ' '.join(filtered_positions), '\n'.join(filtered_genotypes), '',]
    #   Return it!
    return(new_sim)

#####
#   Defining arguments
#####
#   A description of the program
DESCR = """A Python tool to apply a frequency filter to output from an ms
simulation. Reads simulation output from stdin, and outputs in exactly the 
same format, with the filter applied."""
#   Create a new argument parser
Arguments = argparse.ArgumentParser(description=DESCR, add_help=True)
#   Add some arguments
#   This is the panel size
Arguments.add_argument('-p',
                '--panel_size',
                metavar='PANEL_SIZE',
                type=int,
                help='Size of discovery panel',
                required=True)
#   And this is the minimum MAF
Arguments.add_argument('-f',
                '--maf',
                metavar='MAF',
                type=int,
                help='Minimum minor allele count',
                required=True)
#   And parse them
ParsedArgs = Arguments.parse_args()
#   Do some simple argument checking
#   MAF can't be greater than half the size
if ParsedArgs.maf > ParsedArgs.panel_size/2:
    sys.stderr.write('Your MAF limit cannot be greater than half the panel size!\n')
    exit(1)
#   Start a new list to hold simulation data
sim_data = []
lines_read = 0
#   Begin reading through stdin
for lno, line in enumerate(sys.stdin):
    #   The first line contains the command used to run MS
    if lno == 0:
        command = line.strip()
        #   We are interested in the number of chromosomes sampled, which 
        #   determines how many lines the simulation contains
        #   The syntax is
        #   ms nsam nreps [args...]
        nsam = int(command.split(' ')[1])
        #   Then we determine the length of each simulation from this
        #   There is a string for number of segregating sites,
        #   the positions of the segsites
        #   and a newline
        simlength = nsam + 3
    elif lno >= 3:
        #   If the line starts with '//' then we are entering a simulation
        #   Set the counter
        if line.startswith('//'):
            sim_data.append(line.strip())
            lines_read = 1
        if lines_read < simlength:
            sim_data.append(line.strip())
            lines_read += 1
        elif lines_read == simlength:
            #   Append just one more line
            sim_data.append(line.strip())
            #   Then we send the simulation data to be parsed
            new_sim = trim_simulations(sim_data, ParsedArgs.panel_size, ParsedArgs.maf)
            #   And then print it out
            print '\n'.join(new_sim)
            #   Reset the list and the counter
            sim_data = []
            lines_read = 0
