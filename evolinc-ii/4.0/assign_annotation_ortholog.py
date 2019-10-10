__author__ = 'asherkhb' and 'Andrew Nelson'

import sys

#Ensure correct number of arguments are entered.
if len(sys.argv) != 4:
    print("Incorrect usage syntax: Must specify input FASTA, annotation list, and output file name./n")
    print("Use 'rehead.py <Input FASTA> <Annotation List> <Output file name>'")
    exit()

#Establish the inputs.
input_file = sys.argv[1]
annotation_list = sys.argv[2]
output_file = sys.argv[3]

id_dict = {}
result = {}

#Create an output file for writing.
with open(output_file, 'w') as otpt:
    #Open the input FASTA file for processing.
    with open(input_file, 'r') as inpt:
        #Establish annotation list
        annotation_list = open(annotation_list, 'U')
        for line in annotation_list:
            #line_contents = line.strip().split(',')  # Use this line for comma separated
            line_contents = line.strip().split('\t')  # Use this line for tab separated
            id_dict[line_contents[0]] = line_contents[0]
        annotation_list.close()

        #Begin processing FASTA.
        reading = True
        while reading:
            line = inpt.readline()
            #Stop processing at end of input FASTA.
            if line == '':
                reading = False
            #Process FASTA content lines.
            else:
                #Check header lines against list, rename if present in list. Write out appropriate header.
                if line[0] == '>':
                    if "_Known" in line:
                        if "TBH_1" in line:
                            bob,tail = (line.strip('>').strip('\n').strip(' ').split("TBH_1",1))
                            head = '%sTBH_1' % (bob)
                            if head in id_dict:
                                line = '>%s_Homolog%s\n' % (head,tail)
                            else:
                                pass
                        else: 
                            pass
                        otpt.write(line)
                #Write out sequence lines.
                    else:
                        if "TBH_1" in line: 
                            head2 = line.strip('>').strip('\n').strip(' ')
                            if head2 in id_dict:
                                line = '>%s_Homolog\n' % (head2)
                            else:
                                pass
                        else:
                            pass
                        otpt.write(line)
                else:
                    otpt.write(line)
