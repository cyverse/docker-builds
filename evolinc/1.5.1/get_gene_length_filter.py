#!/usr/bin/env python

import sys

infile = sys.argv[1] # input file
outfile = sys.argv[2] # output file

genes = {} # empty list

with open(infile, "rU") as fh_in:
	for line in fh_in:
		line = line.strip()
		if line[0] == ">":
			gene_names = line
			genes[gene_names] = ''
		else:
			genes[gene_names]+=line

for (name,val) in genes.items():
    val = len(val)
    with open(outfile, "a") as fh_out:
        if val > 200:
            fh_out.write(name)
            fh_out.write("\n")
            fh_out.write(genes[name])
            fh_out.write("\n")
