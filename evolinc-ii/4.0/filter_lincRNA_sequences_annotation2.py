#!/usr/bin/env python

import sys

infile = sys.argv[1] # input file for blast output
outfile = sys.argv[2] # output file for filtered blast output

final = list()
with open(infile, 'rU' ) as fh_in:
    for line in fh_in:
        line = line.strip().split()
 	ID = "%s\t%s" % (line[0], line[1])
        pident = line[2]
        evalue = line[10]
        if float(pident) > 99:
            if evalue > '10e-20':
                final.append(ID)

with open(outfile, "w") as fh_out:
    final3 = "\n".join(list(final))
    fh_out.write(final3)
    fh_out.write("\n")
