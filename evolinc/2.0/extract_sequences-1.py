__author__ = 'upendrakumardevisetty'

import sys

accfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

AI_DICT = {}

with open(accfile, "rU") as acc_in:
    for line in acc_in:
        AI_DICT[line[:-1]] = 1

skip = 0

with open(infile, "rU") as fh_in:
    with open(outfile, "w") as fh_out:
        for line in fh_in:
            if line.startswith('>'):
                #line_split = line.split(' ')
                gene = line.strip()
                if gene in AI_DICT:
                    fh_out.write(line)
                    skip = 0
                else:
                    skip = 1
            else:
                if not skip:
                    fh_out.write(line)
