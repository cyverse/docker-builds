#!/usr/bin/env python

import argparse

p = argparse.ArgumentParser()
p.add_argument('-i', '--input', help='Input fasta', required=True)
p.add_argument('-o', '--output', help='Output txt file', required=True)
args = p.parse_args()

genes = {} # empty list

with open(args.input, "rU") as fh_in:
	for line in fh_in:
		line = line.strip()
		if line[0] == ">":
			gene_names = line
			genes[gene_names] = ''
		else:
			genes[gene_names]+=line

for (name,val) in genes.items():
    val = len(val)
    with open(args.output, "a") as fh_out:
       fh_out.write(name[1:])
       fh_out.write("\t")
       fh_out.write(str(val))
       fh_out.write("\n")
