#!/usr/bin/env python
# Upendra Kumar Devisetty
# 11/12/15
# Script to replace the header as well as keep the first id 

import sys

file_in = sys.argv[1]
file_in2 = sys.argv[2]
file_out = sys.argv[3]

result = {}
result2 = {}

with open(file_in) as fh_in:
	with open(file_in2) as fh_in2:
		with open(file_out, 'w') as fh_out: 
			for line in fh_in:
				if line.startswith("CO"):
					line = line.split()
					id = line[1]
					result[id] = ""
				else:
					line = line.split()
					gene = line[1]
					result[id] += gene
					result[id] += " "

			for kee, val in result.items():
				val = val.split()
				if len(val) > 1:
					rename = val[0]
					result2[kee] = rename
				else:
					result2[kee] = val

			for line2 in fh_in2:
				line2 = line2.strip()
				if line2.startswith(">"):
					line2 =line2[1:]
					new = result2[line2]
					new2 = "".join(new)
					fh_out.write(">")
					fh_out.write(str(new2))
					fh_out.write("\n")
				else:
					fh_out.write(line2)
					fh_out.write("\n")
