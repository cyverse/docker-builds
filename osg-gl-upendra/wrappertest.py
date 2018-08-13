#!/usr/bin/env python

import re
import sys

from get_gene_length_filter_3 import get_gene_lengths

if __name__ == "__main__":
    for input_file in sys.argv[1:]:
        output_file = re.sub(r"[.][^.]+$", ".txt", input_file)
        with open(input_file, "rU") as fh_in, open(output_file, "w") as fh_out:
            get_gene_lengths(fh_in, fh_out)
