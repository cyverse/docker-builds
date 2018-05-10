# import scikit-allel
import allel
import sys

# Sys arguments
in_vcf = sys.argv[1]
out_file = sys.argv[2]

callset = allel.read_vcf(in_vcf)
# available keys in vcf file
sorted(callset.keys())
# to get reference from vcf file
callset['variants/REF']
# to get genotype/allel form vcf file
callset['calldata/GT']
# to get genotype infomations in array form
gt = allel.GenotypeArray(callset['calldata/GT'])
# write the output to a file
with open(out_file, 'w') as fh_out:
	fh_out.write(str(gt))