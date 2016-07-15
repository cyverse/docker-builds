__author__ = 'asherkhb'

from subprocess import call
import sys

#Ensure correct number of arguments are entered
if len(sys.argv) != 2:
    print("Incorrect usage syntax: Must specify input list.\nUse 'startup_script.py <input> <query_species>'")
    exit()

#Assign input file from first command argument.
inputfile = sys.argv[1]

with open(inputfile, 'r') as inpt:
    for line in inpt:
        if line == '\n':
            pass
        else:
            line_split = line.split('\t')

            genome = line_split[0]
            sequences = line_split[1]
            gff = line_split[2]
            query_gff = line_split[3]
            query_species = line_split[4]
            query_genome = line_split[5]
            query = "bash /Reciprocal_BLAST.sh -g %s -s %s -f %s -a %s -b %s -c %s" % \
                            (genome, sequences, gff, query_gff, query_species, query_genome)
            if gff == query_gff:
                pass
            else:
                call(query, shell=True)
