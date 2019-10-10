#!/usr/bin/env python
__author__ = 'asherkhb'
# Merge Close Hits
# Version 2.0.2
#
# From an input TSV (of BLAST hits, from evolinc pipeline)
# Label top blast hit, sort by location, then by ID
# Merge close blast hits.
#
# Returns:
#    - sample-instance-report : List of ID's and number of times they occur
#    - <output> : Merged input file.

from sys import argv
from subprocess import call
from os import remove


def renameify(arrayed_line, rep):
    old_name = arrayed_line[0].replace('_TBH', '')
    new_name = "%s_%i" % (arrayed_line[0], rep)
    arrayed_line[0] = new_name
    arrayed_line[8] = arrayed_line[8].replace(old_name, new_name)
    newline = "\t".join(arrayed_line)
    return newline


#Ensure correct number of arguments are entered
if len(argv) != 3:
    print("Incorrect usage syntax: Must specify input and output files.\nUse 'merge_close_hits.py <input> <output>'")
    exit()

#Establish Global Variables
#DONT CHANGE THESE
input_file = argv[1]
input_file_sorted = 'merge_close_hits_input_sorted.txt'
output_file = argv[2]
line_memory = ''
line_memory_array = ["null"]
names = {}


#Label top blast hit, then sort file.
with open(input_file, 'r') as inpt:
    with open(input_file_sorted, 'w') as otpt:
    #Establish a list of IDs
        ids = []
        for line in inpt:
            #Split line, obtain ID
            line_array = line.split('\t')
            line_id = line_array[0]

            #Apply '_TBH' to first instance of
            if line_id not in ids:
                new_line_array = []
                ids.append(line_id)
                for entry in line_array:
                    entry = entry.strip('\n')
                    new_line_array.append(entry)
                new_line_array.append('_TBH')
                new_line = '\t'.join(new_line_array) + '\n'
                otpt.write(new_line)

            else:
                otpt.write(line)

sorting_command = 'sort -k41 --output="%s" %s' % (input_file_sorted, input_file_sorted)
call(sorting_command, shell=True)
echo_command = 'echo >>"%s" %s' % (input_file_sorted, input_file_sorted)
call(echo_command, shell=True)


#Combine Multiple Hits
#Open input file and temporary output file
with open(input_file_sorted, 'r') as inpt:
    with open("scriptpy-temp-out.txt", 'w') as otpt:
        for line in inpt:
            #Skip empty lines
            if line == '\n':
                pass
            #Process data-containing lines
            else:
                #Split line (of TSVs) into an array
                line_array = line.split('\t')

                #Build new line array with newline characters stripped
                new_line_array = []
                for item in line_array:
                    item = item.strip('\n')
                    new_line_array.append(item)
                line_array = new_line_array

                #If two hits are from the same ID, process.
                if line_array[0] == line_memory_array[0] and line_array[1] == line_memory_array[1]:
                    #Establish length of HSP
                    start = int(line_memory_array[3])
                    stop = int(line_array[4])
                    hsp_length = stop - start
                    #Combine HSPs (if appropriate)
                    if int(line_array[9]) + 1000 >= hsp_length > 0:
                        line_memory_array[4] = line_array[4]

                        #Add top blast hit to memory array, if appropriate
                        try:
                            if line_array[10] == '_TBH':
                                line_memory_array.append(line_array[10])
                            else:
                                pass
                        except IndexError:
                            pass

                    #If not appropriate to combine, write out line to temporary output file.
                    else:
                        #Append TBH to ID, if appropriate.
                        try:
                            line_memory_array[0] += line_memory_array[10]
                            del line_memory_array[10]
                        except IndexError:
                            pass

                        #Remove length variable
                        del line_memory_array[9]

                        #Join array into printable line, write to output.
                        printline = "\t".join(line_memory_array) + '\n'
                        otpt.write(printline)

                        #Reset line memory array to current line array.
                        line_memory_array = line_array

                #If two hits are not from same ID, write out line to temporary output file.
                else:
                    if line_memory_array[0] != 'null':
                        #Append TBH to ID, if appropriate.
                        try:
                            line_memory_array[0] += line_memory_array[10]
                            del line_memory_array[10]
                        except IndexError:
                            pass

                        #Remove length variable
                        del line_memory_array[9]

                        #Join array into printable line, write to output.
                        printline = "\t".join(line_memory_array) + '\n'
                        otpt.write(printline)

                    #Reset line memory array to current line array.
                    line_memory_array = line_array


#Rename Duplicates
with open('scriptpy-temp-out.txt', 'r') as inpt:
    with open(output_file, 'w') as otpt:
        for line in inpt:
            #Skip lines without content
            if line == "\n":
                pass

            #Process data-containing lines
            else:
                #Split line (of TSVs) into an array
                line_array = line.split("\t")

                #Check to see if the ID already exists in dictionary
                #If so, increase count and change ID to represent count. Write out line.
                if line_array[0] in names:
                    names[line_array[0]] += 1
                    new_line = renameify(line_array, names[line_array[0]])
                    otpt.write(new_line)

                #If ID does not exist in dictionary, add and set count to 1. Write out line
                else:
                    names[line_array[0]] = 1
                    new_line = renameify(line_array, 1)
                    otpt.write(new_line)

#Remove the temporary output file
remove("scriptpy-temp-out.txt")
remove(input_file_sorted)

#Generate a report of IDs and and number of HSPs
with open('sample-instance-report.txt', 'w') as instance_report:
    for key in names:
        instance = '%i\t%s\n' % (names[key], key)
        instance_report.write(instance)

#Sort the instance report.
call('sort -k2 sample-instance-report.txt --output="sample-instance-report.txt"', shell=True)
