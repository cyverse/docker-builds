__author__ = 'asherkhb'

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
                    head = line.strip('>').strip('\n').strip(' ')
                    if head in id_dict:
                        line = '>%s_Known_lincRNA\n' % (head)
                    else:
                        pass
                    otpt.write(line)
                #Write out sequence lines.
                else:
                    otpt.write(line)
