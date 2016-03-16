#!/usr/bin/env python
import logging
import fnmatch
import os
import numpy as np
from glob import glob
import argparse
import pysam


def get_seq_len_from_bam(samfile):
    temp = []
    for i, dic in enumerate(samfile.header['SQ']):
        try:    
            if args.v:
                print("index: {0} key : {1} value: {2}".format(i,dic['SN'],dic['LN']))
            temp.append({'Seq':dic['SN'], 'Length':dic['LN']})
        except Exception:
            print("i: {0} d: {1}".format(i,dic))
            #print l
            raise
    return temp


def coverage_vectors(contigs_size):
    coverage = {}
    for i in contigs_size:
        coverage[i["Seq"]] = np.zeros(i["Length"])
    return coverage


def parse(samfile,coverage):
    for l in samfile.fetch():
        try:
            if args.pv=="0.7.5":
                # pysam 0.7.5 (for BamM compatibility on OSC)
                begin = l.pos+1
                end = l.pos+l.alen
                coverage[samfile.getrname(l.tid)][begin:end] = 1
            else:
                # pysam 0.8.4
                begin = l.reference_start+1
                end = l.reference_start+l.reference_length
                coverage[samfile.getrname(l.reference_id)][begin:end] = 1
        except Exception:
            print("line: {0}".format(l))
            raise
    return coverage


def write_output(samfile,prop,th,bamout):
    for l in samfile.fetch():
        try:
            if args.pv=="0.7.5":
                # pysam 0.7.5 (for BamM compatibility on OSC)
                prop_g=prop[samfile.getrname(l.tid)]
            else:
                # pysam 0.8.4
                prop_g=prop[samfile.getrname(l.reference_id)]
            if prop_g>th:
                bamout.write(l)
            else:
                if args.v:
                    if args.pv=="0.7.5":
                        # pysam 0.7.5 (for BamM compatibility on OSC)
                        print("we remove read {0} because it matches {1} only covered at {2}% when we ask for at least {3}%".format(l.qname,samfile.getrname(l.tid),prop_g,th))
                    else:
                        # pysam 0.8.4
                        print("we remove read {0} because it matches {1} only covered at {2}% when we ask for at least {3}%".format(l.query_name,samfile.getrname(l.reference_id),prop_g,th))
        except Exception:
            print("line: {0}".format(l))
            raise


def main():
    parser = argparse.ArgumentParser(description='Filter a bam file to keep only reads mapping to sequences covered on at least XX % of their length.')
    # No need of the fasta file anymore, we read sequence length from bam file header. Keep it here in case we have to work with sam files without headers one day
    #parser.add_argument('--fasta','-f', dest='fasta_file', required=True, 
                   #help='Fasta file of the reference sequences')
    parser.add_argument('--bam','-b', dest='bam_file', required=True, 
                   help='input bam file (sorted and indexed)')
    parser.add_argument('--out','-o', dest='output_bam_file', required=True, 
                   help='output bam file')
    parser.add_argument('--th','-t', dest='th', type=int, choices=range(0,100), required=True, 
                   help='minimum coverage of reference sequence (0-100)', metavar="[0-100]")
    parser.add_argument("--verbose",'-v', dest= 'v', help="verbose mode",
                    action="store_true")
    parser.add_argument("--pysam_version",'-pv', dest= 'pv', help="used to specify the version of pysam available (specific variable names for version 0.7.5)",
                    default="0.8.4")
    global args 
    args = parser.parse_args()

    if args.pv=="0.7.5":
        samfile = pysam.Samfile(args.bam_file, "rb") # for older pysam version on OSC
    else:
        samfile = pysam.AlignmentFile(args.bam_file, "rb")

    print("Reading the length of sequences from bam file header")
    contigs_size=get_seq_len_from_bam(samfile)
    print("Getting the coverage tables ready")
    coverage = coverage_vectors(contigs_size)
    print("Getting coverage values")
    coverage = parse(samfile,coverage)
    print("Calculating sequence coverage ratios")
    coverage_prop = {}
    for contig,vector in coverage.items():
        coverage_prop[contig] = np.sum(vector)/float(len(vector))*100
    print("Now writing the output")
    if args.pv=="0.7.5":
        bamout = pysam.Samfile(args.output_bam_file,"wb", template=samfile)
    else:
        bamout = pysam.AlignmentFile(args.output_bam_file,"wb", template=samfile)
    write_output(samfile,coverage_prop,args.th,bamout)
    bamout.close()
    samfile.close()


if __name__ == "__main__":
    output = main()
