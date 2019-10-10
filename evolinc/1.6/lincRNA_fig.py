#!/usr/bin/python

import sys
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

file1 = open(sys.argv[1])
file2 = open(sys.argv[2])
file3 = open(sys.argv[3])

all = []

count_1 = 0
count_2 = 0
count_3 = 0

for line in file1:
	line = line.strip()
	if line.startswith(">"):
		count_1 += 1

all.append(count_1)

for line in file2:
        line = line.strip()
        if line.startswith(">"):
                if re.search("CAGE_PLUS", line):
			count_2 += 1

all.append(count_2)

for line in file3:
        line = line.strip()
        if line.startswith(">"):
		if re.search("known", line):
                	count_3 += 1


all.append(count_3)


labels = ["All lincRNA", "CAGE lincRNA", "Known lincRNA"]
sizes = all
colors = ['yellowgreen', 'mediumpurple', 'lightskyblue'] 
explode = [0, 0.1, 0.1]

plt.pie(sizes,              # data
        explode=explode,
        labels=labels,      # slice labels
	    colors=colors,      # array of colours
        autopct='%1.1f%%',  # print the values inside the wedges
        shadow=True,
        startangle=90
        )
plt.axis('equal')

plt.savefig('lincRNA_piechart')