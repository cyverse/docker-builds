import sys

accfile = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

AI_DICT = {}

with open(accfile, "rU") as acc_in:
    for line in acc_in:
        AI_DICT[line[:-1]] = 1

with open(infile, "rU") as fh_in:
    with open(outfile, "w") as fh_out:
        for line in fh_in:
            line = line.strip()
            line_split = line.split(' ')
            gene = line_split[0]
            if not gene in AI_DICT:
                fh_out.write(gene)
                fh_out.write("\n")