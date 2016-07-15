#!/usr/bin/env python

from os import makedirs, path, listdir
from shutil import move
from sys import argv


def process_file(fasta):
    """Process a FASTA File

    Determines what species represented and the number of times each is represented.

    :param fasta: FASTA file to process.

    :return generate_fasta_entry: Summary of FASTA contents, in TSV format.
    """
    #Establish dictionary of species.
    fasta_spp = {}
    #Open the file for processing.
    file_path = "%s%s" % (directory, fasta)

    #Determine species in file.
    with open(file_path, 'r') as current_file:
        reading = True
        while reading:
            line = current_file.readline()
            if line == '':
                reading = False
            else:
                #Identify header lines.
                if line[0] == '>':
                    #Obtain species from header.
                    line_components = line.split('_')
                    spp = line_components[0].strip('>')
                    #Add species to species list (if not already present).
                    if not spp in fasta_spp:
                        fasta_spp[spp] = 1
                    else:
                        fasta_spp[spp] += 1

    #Move file to folder of most divergent species present.
    moved = 0
    for sps in species_list_r:
        if not moved:
            try:
                if fasta_spp[sps] > 0:
                    new_path = './%s' % sps
                    move(file_path, new_path)
                    moved = 1
            #Skip instances where a possible species is not represented in the FASTA.
            except KeyError:
                pass

    #Return summary file entry.
    return generate_fasta_entry(fasta, species_list_f, fasta_spp)


def generate_fasta_entry(filename, spp_list, fastafile_species_dict):
    """

    :param filename: Name of file being summarized.
    :param spp_list: List of all possible species in FASTA.
    :param fastafile_species_dict: Dictionary of species and occurrences in FASTA.

    :return entry: Summary of FASTA contents, in TSV format.
    """
    #Entry = File Name \t Spp1 Score \t Spp2 Score \t ... SppN Score \n
    entry = filename
    for sps in spp_list:
        if sps in fastafile_species_dict:
            entry += '\t%s' % str(fastafile_species_dict[sps])
        else:
            entry += '\t0'
    entry += '\n'

    return entry


###Begin Functional Script###

#Check for correct number of arguments.
if len(argv) != 2:
    print("Incorrect usage syntax: Must specify species list.")
    print("Use family_division_and_summary.py <species-list>")
    exit()

species_list = argv[1]
species_list_r = []

#Establish FASTA-containing directory, build list of contents.
directory = './'  # Use trailing forward slash.
file_list = listdir(directory)

#Establish species list
with open(species_list, 'r') as spps:
    #Split species file into species list
    species_list_f = spps.read().split('\n')
    #Remove any empty items from list
    if '' in species_list_f:
        species_list_f.remove('')

    #Calculate total species number
    species_number = len(species_list_f)

    #Generate a reversed species list
    for i in range(0, species_number):
        n = (species_number - 1) - i
        species_list_r.append(species_list_f[n])

#Create unique folder for each species.
for sp in species_list_f:
    if not path.exists(sp):
        makedirs(sp)

#Create a summary file and write summarized results in tab-separated form.
with open('summary.txt', 'w') as otpt:
    #Write summary header.
    title = 'File Name'
    for s in species_list_f:
        title += '\t%s' % s
    title += '\n'
    otpt.write(title)

    #Process each file in data directory and write summary of contents to summary.txt.
    for fasta_file in file_list:
        otpt.write(process_file(fasta_file))
