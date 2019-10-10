#!/usr/bin/env python

import sys

def get_gene_lengths(fh_in, fh_out):
    genes = {}

    gene_name = ""
    for line in fh_in:
        line = line.strip()
        if line[0] == ">":
            gene_name = line[1:]
            genes[gene_name] = 0
        elif gene_name != "":
            genes[gene_name] += len(line)

    for (name,val) in genes.items():
        print >>fh_out, "{0}\t{1}".format(name, val)
