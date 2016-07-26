#!/usr/bin/env python
from subprocess import call
import sys
from Bio import SeqIO
import re

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
                
        with open("Homology_Search/List_of_non_identified_query_lincRNAs.txt") as fh_in:
            with open(lincRNA) as fh_in2:    
                for line in fh_in:
                    grab = line.strip()
                    seqiter = SeqIO.parse(fh_in2, 'fasta')
                    new = subject_spp + "." + "Non_identified_query_lincRNAs.fasta"
                    new2 = "Homology_Search" + "/" + new
                    SeqIO.write((seq for seq in seqiter if seq.id in grab), new2, "fasta")

        new = subject_spp + "." + "Non_identified_query_lincRNAs.fasta"
        new2 = "Homology_Search" + "/" + new
        new1 = subject_spp + "." + "Non_identified_query_lincRNAs_singleline.fasta"
        new3 = "Homology_Search" + "/" + new1
        
                                
        with open(new2,"r") as newFile:
            with open(new3,"w") as newFasta:
                sequences = newFile.read()
                sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
                del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
                         # Del removes this empty element.
                for fasta in sequences:
                    try:
                        header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
                    except ValueError:
                        print fasta
                    header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
                    sequence = sequence.replace("\n","") + "\n" # Replace newlines in sequence, remember to add one to the end.
                    newFasta.write(header + sequence)           