#!/usr/bin/env python

import sys
import re

file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as fh_in:
    with open(file_out, 'w') as fh_out:
        for line in fh_in:
		scaffold, program, type, lend, rend, value, strand, empty, morevalues = line.strip().split('\t')
                newlend = int(lend) - 200
                newrend = lend
                if newlend < 0:
                    next
                fh_out.write(scaffold)
                fh_out.write("\t")
                fh_out.write(str(program))
                fh_out.write("\t")
                fh_out.write(str(type))
                fh_out.write("\t")
                fh_out.write(str(newlend))
                fh_out.write("\t")
                fh_out.write(str(newrend))
                fh_out.write("\t")
                fh_out.write(str(value))
                fh_out.write("\t")
                fh_out.write(str(strand))
                fh_out.write("\t")
                fh_out.write(str(empty))
                fh_out.write("\t")
                fh_out.write(str(morevalues))
                fh_out.write("\n")
