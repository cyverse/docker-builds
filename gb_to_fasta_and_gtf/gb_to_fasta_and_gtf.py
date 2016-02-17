#!/usr/bin/python

# Convert a genbank format file of a genome sequence to fasta and
# gtf format. The output files are compatible with cufflinks
# This script was developed for running inside a Docker image and
# is part of a tutorial on RNA-Seq analysis using Cyverse Discovery
# Environment.
#
# usage:
#
# python gb_to_fasta_and_gtf.py <infile.gb> <features.gtf> <sequence.fasta>

import csv
import os
import pdb
from sys import argv
import subprocess

out_gff = open(argv[2], 'w')

seqret_command = ['seqret','-sequence',argv[1],'-outseq',argv[3],'-feature']

subprocess.check_call(seqret_command)

gff_file = os.path.splitext(argv[3])[0] + '.gff'

for row in csv.reader(open(gff_file), delimiter='\t'):
    if row[0].startswith('#'):
        continue
    data = row[8].split(';')
    if row[2] == 'CDS' or row[2] == 'exon':
        geneid = data[1].replace('locus_tag=', '')
        row[8] = 'gene_id "' + geneid + '"; transcript_id "' + geneid+ '";'

        out_gff.write('\t'.join(row) + '\n')

        row[2] = 'exon'
        out_gff.write('\t'.join(row) + '\n')