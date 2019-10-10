#!/usr/bin/env python
# Upendra Kumar Devisetty
# 12/15/15
# Script to parse Bedtools closest ouput

import sys

file = open(sys.argv[1])

file_2 = open(sys.argv[2])

file_3 = open(sys.argv[3], "w")

result = []

dic = {}

for line in file:
    line = line.strip().split("\t")
    if int(line[20]) <= 100 and int(line[20]) > -100:
	test = line[9]
	new = str.split(test, " ")
	final = new[3]
	final2 = final[:-2]
	final3 = final2[1:]
        result.append(final3)


for line in file_2:
    line = line.strip()
    if line.startswith(">"):
        line = line[1:]
        id = line
        dic[id] = " "
    else:
        dic[id] += line

for ele in set(result):
	for kee, val in dic.items():
		new = "_CAGE_PLUS"
		if ele in kee:
			kee2 = kee+new
			dic[kee2] = dic[kee]
			del dic[kee]


for kee2, val2 in dic.items():
	kee2 = ">"+kee2
	val2 = val2.lstrip()
	file_3.write(kee2)
	file_3.write("\n")
	file_3.write(val2)
	file_3.write("\n")
