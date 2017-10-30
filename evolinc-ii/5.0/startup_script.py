#!/usr/bin/env python
from subprocess import call
import sys
from Bio import SeqIO
import re
import getopt


try:
    opts, args = getopt.getopt(sys.argv[1:], 'hb:i:v:', ['inputfile=', 'inputfolder=', 'e_value=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-b', '--inputfile'):
        inputfile_name = arg
    elif opt in ('-i', '--inputfolder'):
        inputfolder_name = arg
    elif opt in ('-v', '--value'):
        value = arg
    else:
        usage()
        sys.exit(2)


with open(inputfile_name, 'r') as inpt:
    for line in inpt:
        if line == '\n':
            pass
        else:
            line = line.strip()
            line_split = line.split('\t')

            subject = line_split[0]
            subject = inputfolder_name + "/" + subject

            lincRNA = line_split[1]
            lincRNA = inputfolder_name + "/" + lincRNA

            query_spp = line_split[2]

            subject_spp = line_split[3]

            try:
                try:
                    subject_gff = line_split[4]
                    subject_gff = inputfolder_name + "/" + subject_gff

                    known_lincRNAs = line_split[5]
                    known_lincRNAs = inputfolder_name + "/" + known_lincRNAs

                    if subject_gff != '':
                        query = "bash /Building_Families.sh -g %s -l %s -q %s -s %s -f %s -k %s -e %s" % \
                                (subject, lincRNA, query_spp, subject_spp, subject_gff, known_lincRNAs, value)

                    else:
                        query = 'bash /Building_Families.sh -g %s -l %s -q %s -s %s -k %s -e %s' % \
                                (subject, lincRNA, query_spp, subject_spp, known_lincRNAs, value)

                except IndexError:
                    subject_gff = line_split[4]
                    subject_gff = inputfolder_name + "/" + subject_gff

                    query = 'bash /Building_Families.sh -g %s -l %s -q %s -s %s -f %s -e %s' % \
                            (subject, lincRNA, query_spp, subject_spp, subject_gff, value)
            except IndexError:
                query = 'bash /Building_Families.sh -g %s -l %s -q %s -s %s -e %s' % \
                        (subject, lincRNA, query_spp, subject_spp, value)
     
        call(query, shell=True)
  
