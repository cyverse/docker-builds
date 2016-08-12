__author__ = 'upendrakumardevisetty'

import sys

infile = sys.argv[1] # input file for blast output
outfile = sys.argv[2] # output file for filtered blast output

final = list()
with open(infile, 'rU' ) as fh_in:
    for line in fh_in:
        line = line.strip()
        line = line.split()
	ID = "%s\t%s" % (line[0], line[1])
        pident = line[2]
        evalue = line[10]
        if float(pident) > 99:
            if evalue > '10e-20':
                final.append(line[0])
                # print line

with open(outfile, "w") as fh_out:
    final_uniq = set(final)
    final2 = "\n".join(list(final_uniq))
    fh_out.write(final2)
