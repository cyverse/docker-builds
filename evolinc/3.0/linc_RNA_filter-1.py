import sys

infile = sys.argv[1] # input file for blast output
outfile = sys.argv[2] # output file for filtered blast output

final = list()
final1 = list()

with open(infile, 'rU') as fh_in:
    for line in fh_in:
        line = line.strip()
        line = line.split()
        if line[0] != line[1]:
            if float(line[2]) >= 99.00:
                final.append(line[1])
        for i in final:
            if i == line[0]:
                final1.append(line[0])


with open(outfile, 'w') as fh_out:
    final_uni = set(final1)
    final2 = "\n".join(list(final_uni))
    fh_out.write(final2)

