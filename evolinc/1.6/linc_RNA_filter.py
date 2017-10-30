import sys

infile = sys.argv[1] # input file for blast output
outfile = sys.argv[2] # output file for filtered blast output

final = list()
with open(infile, 'rU' ) as fh_in:
    for line in fh_in:
        line = line.strip()
        line = line.split()
        if float(line[2]) > 99.00:
            final.append(line[0])
                # print line

with open(outfile, "w") as fh_out:
    final_uniq = set(final)
    final2 = "\n".join(list(final_uniq))
    fh_out.write(final2)
