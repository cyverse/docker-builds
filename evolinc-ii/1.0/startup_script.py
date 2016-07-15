#!/usr/bin/env python
from subprocess import call
import sys

#Ensure correct number of arguments are entered
if len(sys.argv) != 3:
    print("Incorrect usage syntax: Must specify input list.\nUse 'startup_script.py <input>'")
    exit()

#Assign input file from first command argument.
inputfile = sys.argv[1]
inputfolder = sys.argv[2]

with open(inputfile, 'r') as inpt:
    for line in inpt:
        if line == '\n':
            pass
        else:
            line_split = line.split('\t')

            subject = line_split[0]
            subject = inputfolder + "/" + subject

            lincRNA = line_split[1]
            lincRNA = inputfolder + "/" + lincRNA

            query_spp = line_split[2]

            subject_spp = line_split[3]

            try:
                try:
                    subject_gff = line_split[4]
                    subject_gff = inputfolder + "/" + subject_gff

                    known_lincRNAs = line_split[5]
                    known_lincRNAs = inputfolder + "/" + known_lincRNAs

                    if subject_gff != '':
                        query = "sh /Building_Families.sh -g %s -l %s -q %s -s %s -e %s -k %s" % \
                                (subject, lincRNA, query_spp, subject_spp, subject_gff, known_lincRNAs)
                    else:
                        query = 'sh /Building_Families.sh -g %s -l %s -q %s -s %s -k %s' % \
                                (subject, lincRNA, query_spp, subject_spp, known_lincRNAs)

                except IndexError:
                    subject_gff = line_split[4]
                    subject_gff = inputfolder + "/" + subject_gff

                    query = 'sh /Building_Families.sh -g %s -l %s -q %s -s %s -e %s' % \
                            (subject, lincRNA, query_spp, subject_spp, subject_gff)
            except IndexError:
                query = 'sh /Building_Families.sh -g %s -l %s -q %s -s %s' % \
                        (subject, lincRNA, query_spp, subject_spp)
        call(query, shell=True)
