#!/usr/bin/python

import getopt
import sys
import os
import shutil
import collections
from datetime import datetime
import glob
import multiprocessing
from itertools import chain, islice
import numpy as np
import math
from termcolor import colored
import StatisticPair
import common_functions
import gzip
import subprocess


class FilterPair:
    hash_qual_1 = {}
    hash_qual_2 = {}
    hash_qual_1a = {}
    hash_qual_2a = {}
    hash_len_1 = {}
    hash_len_2 = {}
    hash_gc_1 = {}
    hash_gc_2 = {}
    hash_gc_1a = {}
    hash_gc_2a = {}
    tot_Aa_1 = 0
    tot_Ta_1 = 0
    tot_Ga_1 = 0
    tot_Ca_1 = 0
    tot_Na_1 = 0
    tot_Aa_2 = 0
    tot_Ta_2 = 0
    tot_Ga_2 = 0
    tot_Ca_2 = 0
    tot_Na_2 = 0
    tot_A_1 = 0
    tot_T_1 = 0
    tot_C_1 = 0
    tot_G_1 = 0
    tot_N_1 = 0
    tot_A_2 = 0
    tot_T_2 = 0
    tot_C_2 = 0
    tot_G_2 = 0
    tot_N_2 = 0
    n_cont_read_1 = 0
    n_cont_read_1a = 0
    tot_len_a1 = 0
    tot_gc_a1 = 0
    trim_read_1 = 0
    n_cont_read_2 = 0
    short_read_ct_1 = 0
    TotGC1 = 0
    TotLen1 = 0
    count_read_a = 0
    TotLen2 = 0
    TotGC2 = 0
    count_read = 0
    n_read_ct_1 = 0
    ad_trim_1 = 0
    fail_qual_1 = 0
    max_read_len_1 = 0
    min_read_len_1 = 1e-6
    max_trim_len_1 = 0
    min_trim_len_1 = 1e6
    tot_trim_len_1 = 0
    tot_qual_1a_sum = 0
    tot_qual_1_sum = 0
    tot_qual_1_1 = 0
    tot_qual_1_2 = 0
    tot_qual_1_3 = 0
    tot_qual_1_4 = 0
    trim_qual_ct1 = 0
    n_read_ct_2 = 0
    file_p1 = None
    file_p2 = None
    pathname = None
    pathname2 = None
    file_1_path = None
    file_2_path = None
    Trim = 'False'
    CPU = 2
    n_base = 101
    Adapter = None
    per = 1
    qual_thresh = 20
    min_len = 0
    out_file_1 = None
    out_fmt = 'fastq'
    qual_format = 0
    win_size = 5
    ad_trim_2 = 0
    trim_qual_ct2 = 0
    adapter_list = []
    raw_out_dir_1 = None
    raw_out_dir_2 = None
    short_read_ct_2 = 0
    n_cont_read_2a = 0
    tot_len_a2 = 0
    trim_read_2 = 0
    tot_trim_len_2 = 0
    tot_qual_2a_sum = 0
    max_trim_len_2 = 0
    min_trim_len_2 = 1e6
    fail_qual_2 = 0
    tot_gc_a2 = 0
    read_seq_orig_len_1 = 0
    read_seq_orig_len_2 = 0
    tot_qual_2_1 = 0
    tot_qual_2_2 = 0
    tot_qual_2_3 = 0
    tot_qual_2_4 = 0
    tot_qual_2_sum = 0
    out_file_2 = None
    pipeline_flag = None
    out_folder = None
    _gzip_1 = False
    _gzip_2 = False
    no_vis = False

    def __init__(self):
        self.min_size = 0
        try:
            if len(sys.argv) == 1:
                # self.usage()
                # parser.print_help()
                sys.exit(1)
            option, args = getopt.getopt(sys.argv[1:], "ha:b:c:d:e:f:g:i:j:k:l:m:n:p:q:r:s:t:v:", ["p1=",
                                                                                                   "p2=",
                                                                                                   "qfmt=",
                                                                                                   "msz=",
                                                                                                   "nb=",
                                                                                                   "adp=",
                                                                                                   "per=",
                                                                                                   "qthr=",
                                                                                                   "mqual=",
                                                                                                   "o1=",
                                                                                                   "o2=",
                                                                                                   "ofmt=",
                                                                                                   "trim=",
                                                                                                   "wsz=",
                                                                                                   "cpu=",
                                                                                                   "mlk=",
                                                                                                   "pipeline_flag=",
                                                                                                   "out=",
                                                                                                   "no-vis="])
        except getopt.GetoptError:
            print(colored("\nError: wrong input parameter. check command\n", "red"))
            # self.usage()
            # parser.print_help()
            sys.exit(1)

        for opt, value in option:
            if opt in ("-h", "--help"):
                # self.usage()
                sys.exit()
            elif opt in ("-a", "--p1"):
                self.file_p1 = value
                if os.path.exists(self.file_p1):
                    self.pathname = os.path.dirname(self.file_p1)
                    self.pathname = os.path.abspath(self.pathname)
                    self.file_1_path = self.pathname
                else:
                    print(colored("\nError : The given input file does not exist. Rerun the program by giving correct "
                                  "input file\n", "red"))
                    # self.usage()
                    # parser.print_help()
                    sys.exit(1)
            elif opt in ("-b", "--p2"):
                self.file_p2 = value
                if os.path.exists(self.file_p2):
                    self.pathname2 = os.path.dirname(self.file_p2)
                    self.pathname2 = os.path.abspath(self.pathname2)
                    self.file_2_path = self.pathname2
            elif opt in ("-c", "--qfmt"):
                self.qual_format = int(value)
                if self.qual_format == 0:
                    self.qual_format = None
            elif opt in ("-d", "--msz"):
                self.min_size = value
            elif opt in ("-e", "--nb"):
                self.n_base = int(value)
            elif opt in ("-f", "--adp"):
                self.Adapter = value
                self.adapter_list = self.Adapter.split(',')
                if self.Adapter == 'NULL':
                    self.Adapter = None
            elif opt in ("-g", "--per"):
                self.per = float(value)
                if self.per > 1 or self.per < 0:
                    print(colored('Error: Percent threshold (-g, --per) must be between 0-1\n', "red"))
                    sys.exit(1)
            elif opt in ("-i", "--qthr"):
                self.qual_thresh = int(value)
            elif opt in ("-j", "--mqual"):
                self.MinQual = value
            elif opt in ("-k", "--o1"):
                self.out_file_1 = value
            elif opt in ("-l", "--o2"):
                self.out_file_2 = value
            elif opt in ("-m", "--ofmt"):
                self.out_fmt = value
                if self.out_fmt not in ['fastq', 'fasta']:
                    print(colored('Error: Unknown output file format parameter [fastq|fasta]\n', "red"))
                    sys.exit(1)
            elif opt in ('-n', "--trim"):
                self.Trim = value
                if self.Trim not in ['True', 'False']:
                    print(colored('Error: Unknown trim parameter\n', "red"))
                    sys.exit(1)
                elif self.Trim == 'False':
                    self.Trim = None
            elif opt in ("-p", "--wsz"):
                self.win_size = int(value)
            elif opt in ("-q", "--cpu"):
                self.CPU = int(value)
            elif opt in ("-r", "--mlk"):
                self.min_len = int(value)
            elif opt in ("-s", "--pipeline_flag"):
                self.pipeline_flag = value
            elif opt in ("-t", "--out"):
                self.out_folder = value
            elif opt in ("-v", "--no-vis"):
                self.no_vis = value
                if self.no_vis not in ['True', 'False']:
                    print(colored('Error: Unknown visualization parameter [True|False]\n', "red"))
                    sys.exit(1)
                elif self.no_vis == 'False':
                    self.no_vis = None
            else:
                print(colored('Error: in input parameters\n', "red"))
                # self.usage()
                # parser.print_help()
                sys.exit(1)

        if 'gz' in self.file_p1:
            print("["+str(datetime.now())+"] The fastq file is in gz format and uncompressing it...")
            self._gzip_1 = True
            cmd = ["gunzip", self.file_p1]
            p1 = subprocess.Popen(cmd)
            p1.wait()
            if p1.returncode != 0:
                print(colored("Error: during uncompress file", "red"))
                sys.exit(1)

            '''
            read_file = gzip.GzipFile(self.file_p1, 'rb')
            temp = read_file.read()
            read_file.close()
            #   out_file = file(os.path.splitext(os.path.basename(self.file_p1))[0], 'wb')
            out_file = file(self.file_1_path+'/'+os.path.splitext(os.path.basename(self.file_p1))[0], 'wb')
            out_file.write(temp)
            out_file.close()
            self.file_p1 = self.file_1_path+'/'+os.path.splitext(os.path.basename(self.file_p1))[0]
            #   to delete .fq file at end
            #   self.file_1_path = os.path.abspath(self.file_p1)
            '''

        if 'gz' in self.file_p2:
            self._gzip_2 = True
            cmd = ["gunzip", self.file_p2]
            p1 = subprocess.Popen(cmd)
            p1.wait()
            if p1.returncode != 0:
                print(colored("Error: during uncompress file", "red"))
                sys.exit(1)

            '''
            read_file = gzip.GzipFile(self.file_p2, 'rb')
            temp = read_file.read()
            read_file.close()
            #   out_file = file(os.path.splitext(os.path.basename(self.file_p2))[0], 'wb')
            out_file = file(self.file_2_path+'/'+os.path.splitext(os.path.basename(self.file_p2))[0], 'wb')
            out_file.write(temp)
            out_file.close()
            self.file_p2 = self.file_2_path+'/'+os.path.splitext(os.path.basename(self.file_p2))[0]
            #   to delete .fq file at end
            #   self.file_2_path = os.path.abspath(self.file_p2)
            '''

        if self.Trim and self.win_size is None:
            print(colored("\n\nArgument Error: Provide valid window size\n", "red"))
            # self.usage()
            sys.exit(1)
        if self.qual_format is None:
            if self.file_p1 and self.file_p2 and os.path.isfile(self.file_p1) and os.path.isfile(self.file_p2):
                print("["+str(datetime.now())+"] The fastq quality format is not provided therefore detecting the " \
                                              "fastq variant...")
                self.qual_format = self.detect_fastq_variant()
            else:
                print(colored("\nError: Input file can not found\n", "red"))
                # self.usage()
                # parser.print_help()
                sys.exit(1)
            #   print "#################################################################"
        if self.qual_format == 1:
            print(colored("["+str(datetime.now())+"] The fastq quality format is illumina 1.8+", "red"))
            qfmt_verify = self.detect_fastq_variant()
            if qfmt_verify != self.qual_format:
                print(colored("\nError: Wrong quality format\n", "red"))
                sys.exit(1)
        elif self.qual_format == 2:
            print(colored("["+str(datetime.now())+"] The fastq quality format is illumina 1.3+", "red"))
            qfmt_verify = self.detect_fastq_variant()
            if qfmt_verify != self.qual_format:
                print(colored("\nError: Wrong quality format\n", "red"))
                sys.exit(1)
        elif self.qual_format == 3:
            print(colored("["+str(datetime.now())+"] The fastq quality format is Sanger", "red"))
            qfmt_verify = self.detect_fastq_variant()
            if qfmt_verify != self.qual_format:
                print(colored("\nError: Wrong quality format\n", "red"))
                sys.exit(1)
        else:
            print(colored("\nError: Wrong quality format\n", "red"))
            sys.exit(1)

        self.hash_qual_1 = self.get_hash(2, 43)
        self.hash_qual_2 = self.get_hash(2, 43)
        self.hash_qual_1a = self.get_hash(2, 43)
        self.hash_qual_2a = self.get_hash(2, 43)
        self.hash_len_1 = self.get_hash(10, 101)
        self.hash_len_2 = self.get_hash(10, 101)
        self.hash_gc_1 = self.get_hash(10, 101)
        self.hash_gc_2 = self.get_hash(10, 101)
        self.hash_gc_1a = self.get_hash(10, 101)
        self.hash_gc_2a = self.get_hash(10, 101)

    def filter_pair(self):
        #   raw_dir_1 = self.pathname+'/'+'raw_1'
        #   raw_dir_2 = self.pathname+'/'+'raw_2'
        #   self.raw_out_dir_1 = self.pathname+'/'+'raw_out_1'
        #   self.raw_out_dir_2 = self.pathname+'/'+'raw_out_2'
        if self.pipeline_flag == "yes":
            #   pathname will point to output folder pipeline_out
            self.pathname = self.out_folder
            output_dir = self.pathname+'/filtering_out'
            raw_dir_1 = self.pathname+'/'+'raw_1'
            raw_dir_2 = self.pathname+'/'+'raw_2'
            self.raw_out_dir_1 = self.pathname+'/'+'raw_out_1'
            self.raw_out_dir_2 = self.pathname+'/'+'raw_out_2'
        else:
            output_dir = self.pathname+'/'+os.path.splitext(os.path.basename(self.file_p1))[0]+'_filtering_out'
            raw_dir_1 = self.pathname+'/'+'raw_1'
            raw_dir_2 = self.pathname+'/'+'raw_2'
            self.raw_out_dir_1 = self.pathname+'/'+'raw_out_1'
            self.raw_out_dir_2 = self.pathname+'/'+'raw_out_2'

        if os.path.exists(raw_dir_1):
            shutil.rmtree(raw_dir_1)
        if os.path.exists(raw_dir_2):
            shutil.rmtree(raw_dir_2)
        if os.path.exists(self.raw_out_dir_1):
            shutil.rmtree(self.raw_out_dir_1)
        if os.path.exists(self.raw_out_dir_2):
            shutil.rmtree(self.raw_out_dir_2)
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.makedirs(raw_dir_1)
        os.makedirs(raw_dir_2)
        os.makedirs(self.raw_out_dir_1)
        os.makedirs(self.raw_out_dir_2)
        os.makedirs(output_dir)

        print("["+str(datetime.now())+"]" + " Preparing the data for analysis...")
        #   file_p1 = open(self.pathname+'/'+os.path.basename(self.file_p1), 'rU')
        #   file_p2 = open(self.pathname+'/'+os.path.basename(self.file_p2), 'rU')
        file_p1 = open(self.file_1_path+'/'+os.path.basename(self.file_p1), 'rU')
        file_p2 = open(self.file_1_path+'/'+os.path.basename(self.file_p2), 'rU')
        file_p1_basename = os.path.splitext(os.path.basename(self.file_p1))[0]
        file_p2_basename = os.path.splitext(os.path.basename(self.file_p2))[0]
        line_count_1 = sum(1 for line in file_p1)
        line_count_2 = sum(1 for line in file_p2)
        if line_count_1 != line_count_2:
            print(colored("Error: The number of lines are not same in paired end files: Please check the input "
                          "files\n", "red"))
            sys.exit()
        num_split = int(int(line_count_1)/int(self.CPU))
        
        while num_split % 4 != 0:
            num_split += 1
        
        self.split_file(self.file_p1, num_split, raw_dir_1)
        self.split_file(self.file_p2, num_split, raw_dir_2)
        all_files_1 = glob.glob(raw_dir_1+'/*')
        all_files_1 = sorted(all_files_1)
        all_files_2 = glob.glob(raw_dir_2+'/*')
        all_files_2 = sorted(all_files_2)
        Proc = []
        print("["+str(datetime.now())+"]"+" Started Filtering of the reads data...")
        lock = multiprocessing.Lock()
        for i, (file1, file2) in enumerate(zip(all_files_1, all_files_2)):
            #   Without lock
            t = multiprocessing.Process(target=self.filter_data, args=(file1, file2, i))
#           With lock
#           t = multiprocessing.Process(target=filter_data, args=(file, i))
            Proc.append(t)
            t.start()

        for t in Proc:
            t.join()

        StatParse = StatisticPair.StatisticPair(raw_dir_2+'/'+'StatTemp.txt')
        StatParse.stat_pair(self.n_base, self.qual_thresh, self.min_size, self.Adapter, self.Trim, self.min_len,
                            self.file_p1, self.file_p2, self.qual_format, output_dir)
        # for creating plots
        # check for visualization
        # this is important where visualization is not supported
        # by providing -no-vis option, visualization will be off and no error will generated
        if self.no_vis is None:
            StatParse.stat_vis(file_p1_basename, file_p2_basename)

        print("["+str(datetime.now())+"]"+" Finished filtering of data successfully...")

        command_log = open(output_dir+'/Command.log', 'w')
        command_log.write("Pair Filter:\n"+str(sys.argv[0:]))
        command_log.close()

        print(colored("["+str(datetime.now())+"] Output saved in "+output_dir+"\n", "green"))

        if self.pipeline_flag == "yes":
            common_functions.success_flag(self.out_folder+"/success.txt")

        os.chdir(output_dir)
        all_files_out_1 = glob.glob(self.raw_out_dir_1+'/*')
        all_files_out_2 = glob.glob(self.raw_out_dir_2+'/*')
        all_files_out_1 = sorted(all_files_out_1)
        all_files_out_2 = sorted(all_files_out_2)
        if self.out_fmt == "fastq":
            filter_out_1 = open(file_p1_basename+'_Clean.fastq', 'w')
            filter_out_2 = open(file_p2_basename+'_Clean.fastq', 'w')
        elif self.out_fmt == "fasta":
            filter_out_1 = open(file_p1_basename + '_Clean.fasta', 'w')
            filter_out_2 = open(file_p2_basename + '_Clean.fasta', 'w')
        else:
            print(colored('Error: Unknown output file format parameter [fastq|fasta]\n', "red"))
            sys.exit(1)
        self.merge_files(all_files_out_1, filter_out_1)
        self.merge_files(all_files_out_2, filter_out_2)
        shutil.rmtree(raw_dir_1)
        shutil.rmtree(raw_dir_2)
        shutil.rmtree(self.raw_out_dir_1)
        shutil.rmtree(self.raw_out_dir_2)
        if self._gzip_1:
            #   os.remove(self.file_1_path+'/'+self.file_p1)
            os.remove(self.file_p1)
        if self._gzip_2:
            #   os.remove(self.file_1_path+'/'+self.file_p2)
            os.remove(self.file_p2)

    def merge_files(self, all_files_out, filter_out):
        for f in all_files_out:
            if not f.endswith('.txt'):
                temp_file = open(f, 'r')
                for line in temp_file:
                    filter_out.write(line)

    def filter_data(self, fastq_file_1, fastq_file_2, i):
        max_read_len_1 = 0; min_read_len_1 = 1000000
        max_read_len_2 = 0; min_read_len_2 = 1000000
        f1 = open(fastq_file_1, 'rU')
        f2 = open(fastq_file_2, 'rU')
#       l.acquire()
        out_file_1 = str(i)
        out_file_2 = str(i)
#       out_file_1 = open(self.pathname+'/'+'raw_out_1'+'/'+out_file_1, 'w')
#       out_file_2 = open(self.pathname+'/'+'raw_out_2'+'/'+out_file_2, 'w')
        for line1, line2 in zip(f1, f2):
            self.count_read += 1
            header_1_1 = line1.rstrip()
            Header1_2 = line2.rstrip()
            if not header_1_1.startswith('@') or not Header1_2.startswith('@'):
                print("Error: Sequences are not in fastq format\n")
                sys.exit()
            read_seq_1 = next(f1).rstrip()
            read_seq_2 = next(f2).rstrip()
            self.read_seq_orig_len_1 = len(read_seq_1)
            self.read_seq_orig_len_2 = len(read_seq_2)
            header_2_1 = next(f1).rstrip()
            header_2_2 = next(f2).rstrip()
            read_qual_1 = next(f1).rstrip()
            read_qual_2 = next(f2).rstrip()
            #   Count before filtering
            self.tot_A_1 += read_seq_1.count('A')
            self.tot_T_1 += read_seq_1.count('T')
            self.tot_G_1 += read_seq_1.count('G')
            self.tot_C_1 += read_seq_1.count('C')
            self.tot_N_1 += read_seq_1.count('N')
            self.tot_A_2 += read_seq_2.count('A')
            self.tot_T_2 += read_seq_2.count('T')
            self.tot_G_2 += read_seq_2.count('G')
            self.tot_C_2 += read_seq_2.count('C')
            self.tot_N_2 += read_seq_2.count('N')
            if read_seq_1.count('N') > 0:
                self.n_cont_read_1 += 1
            if read_seq_2.count('N') > 0:
                self.n_cont_read_2 += 1
            max_read_len_1, min_read_len_1 = self.find_len(read_seq_1, max_read_len_1, min_read_len_1)
            max_read_len_2, min_read_len_2 = self.find_len(read_seq_2, max_read_len_2, min_read_len_2)
            self.TotLen1 += self.read_seq_orig_len_1
            self.TotLen2 += self.read_seq_orig_len_2
            GC1 = read_seq_1.count('G')+read_seq_1.count('C')
            GC2 = read_seq_2.count('G')+read_seq_2.count('C')
            self.TotGC1 += GC1
            self.TotGC2 += GC2
            if self.qual_format in (1, 3):
                QualFact = 33
                self.process_qual(read_qual_1, GC1, QualFact, 1)
                self.process_qual(read_qual_2, GC2, QualFact, 2)
            elif self.qual_format == 2:
                QualFact = 64
                self.process_qual(read_qual_1, GC1, QualFact, 1)
                self.process_qual(read_qual_2, GC2, QualFact, 2)
#           Actual filtering of data starts here
            if len(read_seq_1) >= int(self.min_size) and len(read_seq_2) >= int(self.min_size):
                if float(read_seq_1.count('N')*100)/len(read_seq_1) <= float(self.n_base) and \
                                        float(read_seq_2.count('N')*100)/len(read_seq_2) <= float(self.n_base):
                    if self.Adapter:
                        read_seq_1, read_qual_1 = self.adapter_trim(read_seq_1, read_qual_1, self.per, 1)
                        read_seq_2, read_qual_2 = self.adapter_trim(read_seq_2, read_qual_2, self.per, 2)
                        if self.Trim:
                            read_seq_1, read_qual_1 = self.qual_trim(read_seq_1, read_qual_1, self.win_size, 1)
                            read_seq_2, read_qual_2 = self.qual_trim(read_seq_2, read_qual_2, self.win_size, 2)
                            self.output_filter_data(header_1_1, header_2_1, read_seq_1, read_qual_1, Header1_2, header_2_2,
                                                    read_seq_2, read_qual_2, out_file_1, out_file_2)
                        else:
                            read_seq_1, read_qual_1 = self.qual_filter(read_seq_1, read_qual_1, 1)
                            read_seq_2, read_qual_2 = self.qual_filter(read_seq_2, read_qual_2, 2)
                            self.output_filter_data(header_1_1, header_2_1, read_seq_1, read_qual_1, Header1_2, header_2_2,
                                                    read_seq_2, read_qual_2, out_file_1, out_file_2)
                    else:
                        if self.Trim:
                            read_seq_1, read_qual_1 = self.qual_trim(read_seq_1, read_qual_1, self.win_size, 1)
                            read_seq_2, read_qual_2 = self.qual_trim(read_seq_2, read_qual_2, self.win_size, 2)
                            self.output_filter_data(header_1_1, header_2_1, read_seq_1, read_qual_1, Header1_2, header_2_2,
                                                    read_seq_2, read_qual_2, out_file_1, out_file_2)
                        else:
                            read_seq_1, read_qual_1 = self.qual_filter(read_seq_1, read_qual_1, 1)
                            read_seq_2, read_qual_2 = self.qual_filter(read_seq_2, read_qual_2, 2)
                            self.output_filter_data(header_1_1, header_2_1, read_seq_1, read_qual_1, Header1_2, header_2_2,
                                                    read_seq_2, read_qual_2, out_file_1, out_file_2)
                else:
                    if float(read_seq_1.count('N')*100)/len(read_seq_1) > float(self.n_base):
                        self.n_read_ct_1 += 1
                    if float(read_seq_2.count('N')*100)/len(read_seq_2) > float(self.n_base):
                        self.n_read_ct_2 += 1
            else:
                if len(read_seq_1) < int(self.min_size):
                    self.short_read_ct_1 += 1
                if len(read_seq_2) < int(self.min_size):
                    self.short_read_ct_2 += 1
        f1.close()
        f2.close()

        StatFile = open('StatTemp.txt', 'a')
        self.hash_gc_1 = collections.OrderedDict(sorted(self.hash_gc_1.items()))
        self.hash_gc_2 = collections.OrderedDict(sorted(self.hash_gc_2.items()))
        hash_gc_1_list = self.hash_gc_1.values()
        hash_gc_2_list = self.hash_gc_1.values()
        self.hash_gc_1a = collections.OrderedDict(sorted(self.hash_gc_1a.items()))
        self.hash_gc_2a = collections.OrderedDict(sorted(self.hash_gc_2a.items()))
        hash_gc_1a_list = self.hash_gc_1a.values()
        hash_gc_2a_list = self.hash_gc_2a.values()
        self.hash_qual_1 = collections.OrderedDict(sorted(self.hash_qual_1.items()))
        self.hash_qual_2 = collections.OrderedDict(sorted(self.hash_qual_2.items()))
        hash_qual_1_list = self.hash_qual_1.values()
        hash_qual_2_list = self.hash_qual_2.values()
        self.hash_qual_1a = collections.OrderedDict(sorted(self.hash_qual_1a.items()))
        self.hash_qual_2a = collections.OrderedDict(sorted(self.hash_qual_2a.items()))
        hash_qual_1a_list = self.hash_qual_1a.values()
        hash_qual_2a_list = self.hash_qual_2a.values()

#       6
        StatFile.write(str(self.count_read)+'\t'+str(self.trim_read_1)+'\t'+str(self.trim_read_2)+'\t' +
                       str(self.short_read_ct_1)+'\t'+str(self.short_read_ct_2)+'\t'+'NA\t')
#       10
        StatFile.write(str(self.TotGC1)+'\t'+str(self.TotGC2)+'\t'+str(self.TotLen1)+'\t'+str(self.TotLen2)+'\t')
#       15
        StatFile.write('NA\t'+str(self.tot_gc_a1)+'\t'+str(self.tot_gc_a2)+'\t'+str(self.tot_len_a1)+'\t' +
                       str(self.tot_len_a2)+'\t')
#       19
        StatFile.write(str(self.count_read_a)+'\t'+str(self.n_read_ct_1)+'\t'+str(self.n_read_ct_2)+'\t'+str(self.ad_trim_1)+'\t')
#       25
        StatFile.write(str(self.ad_trim_2)+'\t'+'NA\t'+str(self.fail_qual_1)+'\t'+str(self.fail_qual_2)+'\t' +
                       str(self.tot_A_1)+'\t'+str(self.tot_A_2)+'\t')
#       29
        StatFile.write(str(self.tot_T_1)+'\t'+str(self.tot_T_2)+'\t'+str(self.tot_G_1)+'\t'+str(self.tot_G_2)+'\t')
#       33
        StatFile.write(str(self.tot_C_1)+'\t'+str(self.tot_C_2)+'\t'+str(self.tot_N_1)+'\t'+str(self.tot_N_2)+'\t')
#       37
        StatFile.write(str(self.tot_Aa_1)+'\t'+str(self.tot_Aa_2)+'\t'+str(self.tot_Ta_1)+'\t'+str(self.tot_Ta_2)+'\t')
#       39
        StatFile.write(str(self.tot_Ga_1)+'\t'+str(self.tot_Ga_2)+'\t')
#       43
        StatFile.write(str(self.tot_Ca_1)+'\t'+str(self.tot_Ca_2)+'\t'+str(self.tot_Na_1)+'\t'+str(self.tot_Na_2)+'\t')
#       47
        StatFile.write(str(self.tot_qual_1_1)+'\t'+str(self.tot_qual_1_2)+'\t'+str(self.tot_qual_1_3)+'\t' +
                       str(self.tot_qual_1_4)+'\t')
#       51
        StatFile.write(str(self.tot_qual_2_1)+'\t'+str(self.tot_qual_2_2)+'\t'+str(self.tot_qual_2_3)+'\t' +
                       str(self.tot_qual_2_4)+'\t')
#       55
        StatFile.write(str(max_read_len_1)+'\t'+str(min_read_len_1)+'\t'+str(max_read_len_2)+'\t'+str(min_read_len_2)+'\t')
#       59
        StatFile.write(str(self.max_trim_len_1)+'\t'+str(self.min_trim_len_1)+'\t'+str(self.max_trim_len_2)+'\t' +
                       str(self.min_trim_len_2)+'\t')
#       71
        StatFile.write(str(self.tot_trim_len_1)+'\t'+str(self.tot_trim_len_2)+'\t' +
                       '\t'.join(str(v) for v in hash_gc_1_list)+'\t')
#       81
        StatFile.write('\t'.join(str(v) for v in hash_gc_2_list)+'\t')
#       91
        StatFile.write('\t'.join(str(v) for v in hash_gc_1a_list)+'\t')
#       101
        StatFile.write('\t'.join(str(v) for v in hash_gc_2a_list)+'\t')
#       105
        StatFile.write(str(self.tot_qual_1_sum)+'\t'+str(self.tot_qual_1a_sum)+'\t'+str(self.tot_qual_2_sum)+'\t' +
                       str(self.tot_qual_2a_sum)+'\t')
#       109
        StatFile.write(str(self.n_cont_read_1)+'\t'+str(self.n_cont_read_1a)+'\t'+str(self.n_cont_read_2)+'\t' +
                       str(self.n_cont_read_2a)+'\t')
#       111
        StatFile.write(str(self.trim_qual_ct1)+'\t'+str(self.trim_qual_ct2)+'\t')
#       132
        StatFile.write('\t'.join(str(v) for v in hash_qual_1_list)+'\t')
#       153
        StatFile.write('\t'.join(str(v) for v in hash_qual_1a_list)+'\t')
#       174
        StatFile.write('\t'.join(str(v) for v in hash_qual_2_list)+'\t')
#       195
        StatFile.write('\t'.join(str(v) for v in hash_qual_2a_list)+'\n')

    def process_qual(self, read_qual_l, GCL, qual_fact_l, flag):
        read_qual_l = list(read_qual_l)
        read_qual_num = map(ord, read_qual_l)
        read_qual_num = [ele-qual_fact_l for ele in read_qual_num]
        if flag == 1:
            self.tot_qual_1_sum += sum(read_qual_num)
        else:
            self.tot_qual_2_sum += sum(read_qual_num)
        avg_qual = np.mean(read_qual_num)
        gc_cont = float((GCL*100)/len(read_qual_l))
        qual_1 = [ele <= 10 for ele in read_qual_num]
        qual_2 = [20 <= ele > 10 for ele in read_qual_num]
        qual_3 = [30 <= ele > 20 for ele in read_qual_num]
        qual_4 = [42 <= ele > 30 for ele in read_qual_num]
#       histogram bins
        if avg_qual == 0:
            avg_qual = 0.1
        if flag == 1:
            self.hash_qual_1[int(math.ceil(avg_qual/2))*2] += 1
            if gc_cont != 0:
                self.hash_gc_1[int(math.ceil(gc_cont/10))*10] += 1
            self.tot_qual_1_1 += len(qual_1)
            self.tot_qual_1_2 += len(qual_2)
            self.tot_qual_1_3 += len(qual_3)
            self.tot_qual_1_4 += len(qual_4)
        else:
            self.hash_qual_2[int(math.ceil(avg_qual/2))*2] += 1
            if gc_cont != 0:
                self.hash_gc_2[int(math.ceil(gc_cont/10))*10] += 1
            self.tot_qual_2_1 += len(qual_1)
            self.tot_qual_2_2 += len(qual_2)
            self.tot_qual_2_3 += len(qual_3)
            self.tot_qual_2_4 += len(qual_4)

    def adapter_trim(self, read_seq_l, read_qual_l, per_l, flag):
        for ad in self.adapter_list:
            ad = list(ad)
            read_seq_l = list(read_seq_l)
            len2 = len(read_seq_l)
            len_ad = len(ad)
            read_seq_l1 = list(read_seq_l)
            read_qual_l = list(read_qual_l)
            iter =0; match = 0
            while iter < len2-len_ad+1 and len2 >= len(ad):
                for i in range(len_ad):
                    if read_seq_l[i] == ad[i]:
                        match += 1
                read_seq_l.pop(0)
                if float(match/len(ad)) >= per_l:
                    del read_seq_l1[iter:iter+len_ad]
                    del read_qual_l[iter:iter+len_ad]
                    del read_seq_l[iter:iter+len_ad]
                    iter = -1
                    len2 = len(read_seq_l)
                iter += 1
                match = 0
            read_seq_l = ''.join(read_seq_l1)
            read_qual_l = ''.join(read_qual_l)
        if len(read_seq_l) < self.read_seq_orig_len_1 and flag == 1:
            self.ad_trim_1 += 1
        elif len(read_seq_l) < self.read_seq_orig_len_2 and flag == 2:
            self.ad_trim_2 += 1
        return read_seq_l, read_qual_l

    def qual_trim(self, read_seq_l, read_qual_l, win_size_l, flag):
        read_seq_l = list(read_seq_l)
        read_qual_l = list(read_qual_l)
        lenb = len(read_seq_l)
        read_qual_num_l = map(ord, read_qual_l)
        temp_read_qual_num_l = []
        temp_read_seq_l = []
        temp_read_qual_l = []
        iter = 0
        if self.qual_format == 1 or self.qual_format == 3:
            read_qual_num_l = [ele-33 for ele in read_qual_num_l]
        if self.qual_format == 2:
            read_qual_num_l = [ele-64 for ele in read_qual_num_l]
        while iter < (len(read_qual_num_l) - win_size_l)+1:
            trim = read_qual_num_l[:win_size_l]
            if np.mean(trim) <= int(self.qual_thresh):
                read_qual_num_l = read_qual_num_l[win_size_l:]
                read_seq_l = read_seq_l[win_size_l:]
                read_qual_l = read_qual_l[win_size_l:]
                iter -= 1
            else:
                temp_read_qual_num_l.extend(read_qual_num_l[:win_size_l])
                read_qual_num_l = read_qual_num_l[win_size_l:]
                temp_read_seq_l.extend(read_seq_l[:win_size_l])
                read_seq_l = read_seq_l[win_size_l:]
                temp_read_qual_l.extend(read_qual_l[:win_size_l])
                read_qual_l = read_qual_l[win_size_l:]
                iter -= 1
            iter += 1
        temp_read_qual_num_l.extend(read_qual_num_l)
        temp_read_seq_l.extend(read_seq_l)
        temp_read_qual_l.extend(read_qual_l)
        if len(temp_read_seq_l) > 0 and np.mean(temp_read_qual_num_l) > int(self.qual_thresh):
            if len(temp_read_seq_l) < lenb and flag == 1:
                self.trim_qual_ct1 += 1
            elif len(temp_read_seq_l) < lenb and flag == 2:
                self.trim_qual_ct2 += 1
            return ''.join(temp_read_seq_l), ''.join(temp_read_qual_l)
        else:
            return None, None

    def output_filter_data(self, header_1_1_l, header_2_1_l, read_seq_1_l, read_qual_1_l, header_1_2_l, header_2_2_l, 
                           read_seq_2_l, read_qual_2_l, out_file_1_l, out_file_2_l):
        out_file_1_l = open(self.pathname+'/'+'raw_out_1'+'/'+out_file_1_l, 'a')
        out_file_2_l = open(self.pathname+'/'+'raw_out_2'+'/'+out_file_2_l, 'a')
        if out_file_1_l and out_file_2_l and read_seq_1_l and read_seq_2_l and len(read_seq_1_l) > int(self.min_len) and \
                        len(read_seq_2_l) > int(self.min_len) and self.out_fmt == 'fastq':
            self.count_read_a += 1
            out_file_1_l.write(header_1_1_l+'\n'+read_seq_1_l+'\n'+header_2_1_l+'\n'+read_qual_1_l+'\n')
            out_file_2_l.write(header_1_2_l+'\n'+read_seq_2_l+'\n'+header_2_2_l+'\n'+read_qual_2_l+'\n')
            self.output_filter_data_sub(read_seq_1_l, read_qual_1_l, read_seq_2_l, read_qual_2_l)
        elif out_file_1_l and out_file_2_l and read_seq_1_l and read_seq_2_l and len(read_seq_1_l) > int(self.min_len) and \
                        len(read_seq_2_l) > int(self.min_len) and self.out_fmt == 'fasta':
            self.count_read_a += 1
            out_file_1_l.write('>'+header_1_1_l+'\n'+read_seq_1_l+'\n')
            out_file_2_l.write('>'+header_1_2_l+'\n'+read_seq_2_l+'\n')
            self.output_filter_data_sub(read_seq_1_l, read_qual_1_l, read_seq_2_l, read_qual_2_l)

    def output_filter_data_sub(self, read_seq_1_l, read_qual_1_l, read_seq_2_l, read_qual_2_l):
        self.tot_Aa_1 += read_seq_1_l.count('A')
        self.tot_Ta_1 += read_seq_1_l.count('T')
        self.tot_Ga_1 += read_seq_1_l.count('G')
        self.tot_Ca_1 += read_seq_1_l.count('C')
        self.tot_Na_1 += read_seq_1_l.count('N')
        self.tot_Aa_2 += read_seq_2_l.count('A')
        self.tot_Ta_2 += read_seq_2_l.count('T')
        self.tot_Ga_2 += read_seq_2_l.count('G')
        self.tot_Ca_2 += read_seq_2_l.count('C')
        self.tot_Na_2 += read_seq_2_l.count('N')
        if read_seq_1_l.count('N') > 0:
            self.n_cont_read_1a += 1
        if read_seq_2_l.count('N') > 0:
            self.n_cont_read_2a += 1
        self.tot_len_a1 += len(read_seq_1_l)
        self.tot_len_a2 += len(read_seq_2_l)
        self.tot_gc_a1 += read_seq_1_l.count('G')+read_seq_1_l.count('C')
        self.tot_gc_a2 += read_seq_2_l.count('G')+read_seq_2_l.count('C')
        gc_cont1 = float(((read_seq_1_l.count('G')+read_seq_1_l.count('C'))*100)/len(read_seq_1_l))
        gc_cont2 = float(((read_seq_2_l.count('G')+read_seq_2_l.count('C'))*100)/len(read_seq_2_l))
        if gc_cont1 != 0:
            self.hash_gc_1a[int(math.ceil(gc_cont1/10))*10] += 1
        if gc_cont2 != 0:
            self.hash_gc_2a[int(math.ceil(gc_cont2/10))*10] += 1
        if len(read_seq_1_l) < self.read_seq_orig_len_1:
            self.trim_read_1 += 1
            self.tot_trim_len_1 += len(read_seq_1_l)
        if len(read_seq_2_l) < self.read_seq_orig_len_2:
            self.trim_read_2 += 1
            self.tot_trim_len_2 += len(read_seq_2_l)
        self.max_trim_len_1, self.min_trim_len_1 = self.find_len(read_seq_1_l, self.max_trim_len_1, self.min_trim_len_1)
        self.max_trim_len_2, self.min_trim_len_2 = self.find_len(read_seq_2_l, self.max_trim_len_2, self.min_trim_len_2)
        read_qual_1_l = list(read_qual_1_l)
        read_qual_2_l = list(read_qual_2_l)
        read_qual_num_l1 = map(ord, read_qual_1_l)
        read_qual_num_l2 = map(ord, read_qual_2_l)
        if self.qual_format in (1, 3):
            read_qual_num_l1 = [ele-33 for ele in read_qual_num_l1]
            read_qual_num_l2 = [ele-33 for ele in read_qual_num_l2]
        else:
            read_qual_num_l1 = [ele-64 for ele in read_qual_num_l1]
            read_qual_num_l2 = [ele-64 for ele in read_qual_num_l2]
        avg_qual_1a = np.mean(read_qual_num_l1)
        avg_qual_2a = np.mean(read_qual_num_l2)
        self.hash_qual_1a[int(math.ceil(avg_qual_1a/2))*2] += 1
        self.hash_qual_2a[int(math.ceil(avg_qual_2a/2))*2] += 1
        self.tot_qual_1a_sum += sum(read_qual_num_l1)
        self.tot_qual_2a_sum += sum(read_qual_num_l2)

    def qual_filter(self, read_seq_l, read_qual_l, flag):
        read_qual_l = list(read_qual_l)
        read_qual_num_l = map(ord, read_qual_l)
        if self.qual_format in (1, 3):
            read_qual_num_l = [ele-33 for ele in read_qual_num_l]
            if np.mean(read_qual_num_l) >= int(self.qual_thresh):
                return read_seq_l, ''.join(read_qual_l)
            else:
                if flag == 1:
                    self.fail_qual_1 += 1
                else:
                    self.fail_qual_2 += 1
                return None, None
        elif self.qual_format == 2:
            read_qual_num_l = [ele-64 for ele in read_qual_num_l]
            if np.mean(read_qual_num_l) >= self.qual_thresh:
                return read_seq_l, ''.join(read_qual_l)
            else:
                if flag == 1:
                    self.fail_qual_1 += 1
                else:
                    self.fail_qual_2 += 1
            return None, None

    def detect_fastq_variant(self):
        Count = 0
        check = []
        File = open(self.file_p1, 'rU')
        for line in File:
            id =line.rstrip()
            if not id.startswith('@'):
                print("Error: Sequences are not in fastq format\n")
                sys.exit()
            next(File).rstrip()
            next(File)
            asc = next(File).rstrip()
            asc_list = list(asc)
            asc_list = list(map(ord, asc_list))
            Min = min(asc_list)
            Max = max(asc_list)
            check.append(Min)
            check.append(Max)
            Count += 1
            if Count == 40000:
                break
        File.close()
        Min = min(check)
        Max = max(check)
        if 64 > Min >= 33 and Max == 74:
            return 1
        elif Min >= 64 and 74 < Max <= 104:
            return 2
        elif 64 > Min >= 33 and Max <= 73:
            return 3

    def split_file(self, File, num_lines, Dir):
        #   File = self.pathname+'/'+os.path.basename(File)
        File = self.file_1_path+'/'+os.path.basename(File)
        os.chdir(Dir)
        with open(File) as MainFile:
            for i, lines in enumerate(self.file_parts(MainFile, num_lines)):
                file_split = '{}'.format(i)
                with open(file_split, 'w') as f:
                    f.writelines(lines)
    
    @staticmethod
    def find_len(seq, max_read_len_1_l, min_read_len_1_l):
        if len(seq) >= max_read_len_1_l:
            max_read_len_1_l = len(seq)
        if len(seq) < min_read_len_1_l:
            min_read_len_1_l = len(seq)
        return max_read_len_1_l, min_read_len_1_l

    @staticmethod
    def get_hash(start, end):
        d = {}
        for v in range(start, end, start):
            d[v] = 0
        return d

    @staticmethod
    def file_parts(iterable, n):
        iterable = iter(iterable)
        while True:
            yield chain([next(iterable)], islice(iterable, n-1))

    @staticmethod
    def usage():
        print("\nusage: srap -filter-p [options]")

        print("\nOptions:")
        print("  "+"-a, --p1 INT {:>44}".format("Input file left read (.fastq, .fq)"))
        print("  "+"-b, --p2 INT {:>45}".format("Input file right read (.fastq, .fq)"))
        print("  "+"-d, --msz INT {:>60}".format("Filter the reads which are lesser than minimum size"))
        print("  "+"-c, --qfmt INT {:>28}".format("Quality value format"))
        print("{:>40}".format("1= Illumina 1.8"))
        print("{:>40}".format("2= Illumina 1.3"))
        print("{:>35}".format("3= Sanger\n"))
        print("  "+"-e, --nb INT {:>51}".format("Filter the reads with more N (Value in %)"))
        print("  "+"-f, --adp STRING {:>81}".format("Trim the adapter sequence and truncate the read "
                                                    "sequence [Adapter sequence]"))
        print("  "+"-g, --per FLOAT {:>124}".format("Truncate the read sequence if it matches to adapter sequence "
                                                    "equal or more than given percent (0.0-1.0) [default=0.9]"))
        print("  "+"-i, --qthr INT {:>113}".format("Filter the read sequence if average quality of "
                                                   "bases in reads is lower than threshold (1-40) [default:20]"))
#        print "  "+"-j, --mqual INT {:>122}".format("filter the reads if given percentage of bases having lower "
#                                                    "quality values lower than threshold (1-100) [default:20]")
        print("  "+"-n, --trim BOOLEAN {:>142}".format("If trim option set to true, the reads with low "
                                                       "quality (as defined by option -qthr) will be trimmed instead"
                                                       " of discarding [default: false]"))
        print("  "+"-p, --wsz INT {:>118}".format("The window size for trimming (5\'->3\') the reads. This option"
                                                  " should always set when -trim option is defined. "))
        print("{:>101}".format("The recommended windowsize for best result should be between 2-5 [default:5]"))
        print("  "+"-r, --mlk INT {:>118}".format("The Minimum length of the reads to retain after trimming"))
        print("  "+"-q, --cpu {:>38}".format("Number of CPU [default:2]"))
        print("  "+"-m, --ofmt {:>60}".format("Output file format (fastq/fasta) [default:fastq]"))
        print("  "+"-m, --no-vis BOOLEAN {:>51}".format("No figures will be produced [yes|no] [default:no]"))
        print("  "+"-h, --help {:>37}".format("Print this help message\n\n"))


if __name__ == '__main__':
    Ob = FilterPair()
    try:
        Ob.filter_pair()
    except KeyboardInterrupt:
        pass
        print(colored("\nProgram is terminated by user", "red"))
        sys.exit(1)
