#!/usr/bin/env python

import pandas as pd
import sys

file_in = open(sys.argv[1])
file_in2 = open(sys.argv[2])
file_out = sys.argv[3]

sites = []

for line in file_in:
        line = line.strip()
        if line.startswith(">"):
                line = line[1:]
                gene = str.split(line, ".")
                gene = gene[0]
                sites.append(gene)

pattern = "|". join(sites)

df = pd.read_csv(file_in2, sep = '\t', header = None, engine = 'python')

mask = df[8].str.contains(pattern)

df.loc[mask,1] = 'lincRNA'

df.to_csv(file_out, sep = "\t", index = False, header = False)