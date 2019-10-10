#!/usr/bin/env python

import sys
import pandas as pd

infile1 = sys.argv[1] 
infile2 = sys.argv[2] 
outfile = sys.argv[3] 

df = pd.read_csv(infile2, sep="\t")

final1 = []
final2 = []

with open(infile1) as fh_in1:
	for line in fh_in1:
		line1 = line.strip().split('_')
		line2 = line.strip().split('\t')
		final1.append(line2[0][5:])
		final2.append(line2[1])
	new = line1[0]
	test = df[["id",new]] 	
	dic1 = dict(zip(test.id,test.iloc[:,1]))

dic2 = dict(zip(final1,final2))

with open(outfile, 'w') as fh_out:	
	for key in dic1:
		if key in dic2:
			fh_out.write(key)
			fh_out.write('\t')
			comb = dic1[key] + "_" + dic2[key]
			fh_out.write(comb)
			fh_out.write('\n')			
		else:
			fh_out.write(key)
			fh_out.write('\t')
			fh_out.write(dic1[key])
			fh_out.write('\n')

df2 = pd.read_csv(outfile, sep = '\t', names=["ID", new]) 
df2.to_csv(outfile, sep = '\t', na_rep="NA", index=False)
