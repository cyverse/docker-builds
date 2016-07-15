__author__ = 'upendrakumardevisetty'
#edited by andrew.d.l.nelson@gmail.com

import sys

listfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

AI_DICT = {}

with open(listfile, "rU") as list_in:
    for line in list_in:
        AI_DICT[line[:-1]] = 1

skip = 0

with open(infile, "rU") as fh_in:
    with open(outfile, "w") as fh_out:
        for line in fh_in:
            if line.startswith('>'):
                line_split = line.split('=')
                gene = line_split[0]
                if gene in AI_DICT:
                    fh_out.write(line)
                    skip = 0
                else:
                    skip = 1
            else:
                if not skip:
                    fh_out.write(line)
