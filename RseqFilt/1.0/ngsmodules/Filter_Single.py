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
import StatisticSingle
import common_functions
import subprocess


class FilterSingle:
    hash_qual_1 = {}
    hash_qual_1a = {}
    hash_len_1 = {}
    hash_gc_1 = {}
    hash_gc_1a = {}
    tot_Aa_1 = 0
    tot_Ta_1 = 0
    tot_Ga_1 = 0
    tot_Ca_1 = 0
    tot_Na_1 = 0
    tot_A_1 = 0
    tot_T_1 = 0
    tot_C_1 = 0
    tot_G_1 = 0
    tot_N_1 = 0
    n_cont_read_1 = 0
    n_cont_read_1a = 0
    tot_len_a1 = 0
    tot_gc_a1 = 0
    trim_read_1 = 0
    short_read_ct_1 = 0
    tot_gc_1 = 0
    tot_len_1 = 0
    count_read_a = 0
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
    trim_qual_ct_1 = 0
    n_read_ct_2 = 0
    file_1 = None
    pathname = None
    Trim = 'False'
    CPU = 2
    n_base = 101
    adapter = None
    Per = 1
    qual_thresh = 20
    min_len = 0
    out_file_1 = None
    out_fmt = 'fastq'
    qual_format = 0
    win_size = 5
    ad_trim_2 = 0
    adapter_list = []
    raw_out_dir1 = None
    read_seq_orig_len_1 = 0
    min_size = 0
    low_qual_count = 0
    pipeline_flag = "no"
    out_folder = None
    order = None
    file_1_path = None
    #   for uncompressed gz file
    file_1_path_gz = None
    no_vis = False

    def __init__(self):
        try:
            if len(sys.argv) == 1:
                # self.usage()
                parser.print_help()
                sys.exit(1)
            option, args = getopt.getopt(sys.argv[1:], "ha:b:c:d:e:f:g:i:j:k:l:m:n:p:q:r:s:t:u:v:", ["p1=",
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
                                                                                                     "order=",
                                                                                                     "no-vis="])
        except getopt.GetoptError:
            print(colored("\nError: wrong input parameter. check command\n", "red"))
            # self.usage()
            sys.exit(1)

        for opt, value in option:
            if opt in ('-h', "--help"):
                self.usage()
                sys.exit(1)
            elif opt in ("-a", "--p1"):
                self.file_1 = value
                if os.path.exists(self.file_1):
                    if 'gz' in self.file_1:
                        self.file_1 = os.path.abspath(os.path.splitext(self.file_1)[0])
                        self.pathname = os.path.dirname(self.file_1)
                        self.pathname = os.path.abspath(self.pathname)
                        #   for gz uncompressed file
                        self.file_1_path_gz = self.file_1_path
                    else:
                        self.file_1 = os.path.abspath(self.file_1)
                        self.pathname = os.path.dirname(self.file_1)
                        self.pathname = os.path.abspath(self.pathname)
                        self.file_1_path = self.pathname
                else:
                    print(colored("\nError : The given input file does not exist. Rerun the program by giving correct "
                                  "input file\n", "red"))
                    # self.usage()
                    # parser.print_help()
                    sys.exit(1)
            elif opt in ("-c", "--qfmt"):
                self.qual_format = int(value)
                if self.qual_format == 0:
                    self.qual_format = None
            elif opt in ("-d", "--msz"):
                self.min_size = value
            elif opt in ("-e", "--nb"):
                self.n_base = int(value)
            elif opt in ("-f", "--adp"):
                self.adapter = value
                self.adapter_list = self.adapter.split(',')
                if self.adapter=='NULL':
                    self.adapter = None
            elif opt in ("-g", "--per"):
                self.Per = float(value)
                if self.Per > 1 or self.Per < 0:
                    print(colored('Error: Percent threshold (-g, --per) must be between 0-1\n', "red"))
                    sys.exit(1)
            elif opt in ("-i", "--qthr"):
                self.qual_thresh = int(value)
            elif opt in ('-n', "--trim"):
                self.Trim = value
                if self.Trim not in ['True', 'False']:
                    print(colored('Error: Unknown trim parameter\n', "red"))
                    sys.exit(1)
                elif self.Trim == 'False':
                    self.Trim = None
            elif opt in ("-j", "--mqual"):
                self.min_qual = value
            elif opt in ("-k", "--o1"):
                self.out_file_1 = value
            elif opt in ("-m", "--ofmt"):
                self.out_fmt = value
                if self.out_fmt not in ['fastq', 'fasta']:
                    print(colored('Error: Unknown output file format parameter [fastq|fasta]\n', "red"))
                    sys.exit(1)
            elif opt in ("-p", "--wsz"):
                self.win_size = int(value)
            elif opt in ("-q", "--cpu"):
                self.CPU = value
            elif opt in ("-r", "--mlk"):
                self.min_len = value
            elif opt in ("-r", "--mlk"):
                self.min_len = value
            elif opt in ("-s", "--pipeline_flag"):
                self.pipeline_flag = value
            elif opt in ("-t", "--out"):
                self.out_folder = value
            elif opt in ("-u", "--order"):
                self.order = value
            elif opt in ("-v", "--no-vis"):
                self.no_vis = value
                if self.no_vis not in ['True', 'False']:
                    print(colored('Error: Unknown visualization parameter [True|False]\n', "red"))
                    sys.exit(1)
                elif self.no_vis == 'False':
                    self.no_vis = None
            elif opt in ("-b", "--p2"):
                print(colored("Error: use -filter-p option for paired-end data", "red"))
                sys.exit(1)
            else:
                sys.exit(1)

        if 'gz' in self.file_1:
            print("["+str(datetime.now())+"] The fastq file is in gz format and uncompressing it...")
            self._gzip_1 = True
            cmd = ["gunzip", self.file_1]
            p1 = subprocess.Popen(cmd)
            p1.wait()
            if p1.returncode != 0:
                print(colored("Error: during uncompress file", "red"))
                sys.exit(1)
            #   read_file = gzip.GzipFile(self.file_1, 'rb')
            #   temp = read_file.read()
            #   read_file.close()
            #   out_file = file(os.path.splitext(os.path.basename(self.file_1))[0], 'wb')
            #   out_file = file(self.file_1_path_gz, 'wb')
            #   out_file.write(temp)
            #   out_file.close()
            #   self.file_1 = os.path.splitext(os.path.basename(self.file_1))[0]
            #   self.file_1 = self.file_1_path_gz
            self.file_1 = self.file_1_path
        if self.Trim and self.win_size is None:
            print(colored("\nArgument Error: Provide valid Window size\n", "red"))
            sys.exit(1)
            # self.usage()
        if self.qual_format is None:
            if self.file_1 and os.path.isfile(self.file_1):
                print("["+str(datetime.now())+"] The fastq quality format is not provided therefore detecting the " \
                                              "fastq variant...")
                self.qual_format = self.detect_fastq_variant()
            else:
                print(colored("Error: Input File Can not found", "red"))
                # self.usage()
                sys.exit()

        if self.qual_format == 1:
            print(colored("[" + str(datetime.now()) + "] The fastq quality format is illumina 1.8+", "red"))
            qfmt_verify = self.detect_fastq_variant()
            if qfmt_verify != self.qual_format:
                print(colored("\nError: Wrong quality format\n", "red"))
                sys.exit(1)
        elif self.qual_format == 2:
            print(colored("[" + str(datetime.now()) + "] The fastq quality format is illumina 1.3+", "red"))
            qfmt_verify = self.detect_fastq_variant()
            if qfmt_verify != self.qual_format:
                print(colored("\nError: Wrong quality format\n", "red"))
                sys.exit(1)
        elif self.qual_format == 3:
            print(colored("[" + str(datetime.now()) + "] The fastq quality format is Sanger", "red"))
            qfmt_verify = self.detect_fastq_variant()
            if qfmt_verify != self.qual_format:
                print(colored("\nError: Wrong quality format\n", "red"))
                sys.exit(1)
        else:
            print(colored("\nError: Wrong quality format\n", "red"))
            sys.exit(1)

        self.hash_qual_1 = self.get_hash(2, 43)
        self.hash_qual_1a = self.get_hash(2, 43)
        self.hash_len_1 = self.get_hash(10, 101)
        self.hash_gc_1 = self.get_hash(10, 101)
        self.hash_gc_1a = self.get_hash(10, 101)

    def filter_single(self):
        '''
        raw_dir1 = self.pathname+'/'+'Raw1'
        self.raw_out_dir1 = self.pathname+'/'+'raw_out_1'
        if self.pipeline_flag == "yes":
            #   output_dir = self.pathname+'/'+'Pipeline_Output'+'/filtering_out'
            output_dir = self.out_folder+'/filtering_out'
        else:
            output_dir = self.pathname+'/'+'filtering_out'
        '''

        if self.pipeline_flag == "yes":
            #   pathname will point to output folder pipeline_out
            self.pathname = self.out_folder
            output_dir = self.pathname+'/filtering_out'
            raw_dir1 = self.pathname+'/'+'raw1'
            self.raw_out_dir1 = self.pathname+'/'+'raw_out_1'
        else:
            # output_dir = self.pathname+'/'+'filtering_out'
            output_dir = self.pathname + '/' + os.path.splitext(os.path.basename(self.file_1))[0] + '_filtering_out'
            # self.file_1_path
            raw_dir1 = self.pathname+'/'+'raw1'
            self.raw_out_dir1 = self.pathname+'/'+'raw_out_1'

        if os.path.exists(raw_dir1):
            shutil.rmtree(raw_dir1)
        if os.path.exists(self.raw_out_dir1):
            shutil.rmtree(self.raw_out_dir1)
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.makedirs(raw_dir1)
        os.makedirs(self.raw_out_dir1)
        os.makedirs(output_dir)

        print("["+str(datetime.now())+"] Preparing the data for analysis...")

        if 'gz' in self.file_1:
            file_s = open(self.file_1_path, 'rU')
            file_s_basename = os.path.splitext(os.path.basename(self.file_1_path))[0]
        else:
            file_s = open(self.file_1, 'rU')
            file_s_basename = os.path.splitext(os.path.basename(self.file_1))[0]

        line_count = sum(1 for line in file_s)
        num_split = int(int(line_count)/int(self.CPU))

        while num_split % 4 != 0:
            num_split += 1

        #   Split file into several parts contain num_split lines
        if 'gz' in self.file_1:
            self.split_file(self.file_1_path, num_split, raw_dir1)
        else:
            self.split_file(self.file_1, num_split, raw_dir1)
        #   Read all file names and store in a list
        all_file_s = glob.glob(raw_dir1+'/*')
        all_file_s = sorted(all_file_s)
        Proc = []
        print("["+str(datetime.now())+"] Started Filtering of the reads data...")
        lock = multiprocessing.Lock()
        for i, file in enumerate(all_file_s):
            #   t = multiprocessing.Process(target=filter_data, args=(file, i))
            t = multiprocessing.Process(target=self.filter_data, args=(file, i))
            Proc.append(t)
            t.start()

        for t in Proc:
            t.join()

        print("["+str(datetime.now())+"] Finished filtering of data successfully ...")
        print(colored("["+str(datetime.now())+"] Output saved in "+output_dir+"\n", "green"))

        if self.pipeline_flag == "yes":
            common_functions.success_flag(self.out_folder+"/success.txt")

        os.chdir(output_dir)
        all_file_s_out = glob.glob(self.raw_out_dir1+'/*')
        all_file_s_out = sorted(all_file_s_out)
        if self.out_fmt == "fastq":
            filter_out = open(file_s_basename+'_Clean.fastq', 'w')
        elif self.out_fmt == "fasta":
            filter_out = open(file_s_basename + '_Clean.fasta', 'w')
        else:
            print(colored('Error: Unknown output file format parameter [fastq|fasta]\n', "red"))
            sys.exit(1)

        for f in all_file_s_out:
            if not f.endswith('.txt'):
                temp_file = open(f, 'r')
                for line in temp_file:
                    filter_out.write(line)

        stat_parse = StatisticSingle.StatisticSingle(self.raw_out_dir1+'/'+'StatTemp.txt')
        stat_parse.stat_single(self.n_base, self.qual_thresh, self.min_size, self.adapter, self.Trim, self.min_len,
                              self.file_1, self.qual_format, output_dir)

        # check for visualization
        # this is important where visualization is not supported
        # by providing -no-vis option, visualization will be off and no error will generated
        if self.no_vis is None:
            stat_parse.stat_vis(file_s_basename)

        shutil.rmtree(raw_dir1)
        shutil.rmtree(self.raw_out_dir1, ignore_errors=True)

    def split_file(self, File, NumLines, Dir):
        #   File = self.pathname+'/'+os.path.basename(File)
        #   change directory
        os.chdir(Dir)
        with open(File) as MainFile:
            for i, lines in enumerate(self.file_parts(MainFile, NumLines)):
                file_split = '{}'.format(i)
                with open(file_split, 'w') as f:
                    f.writelines(lines)

    def filter_data(self, fastq_file, i):
        max_read_len_1 = 0; min_read_len_1 = 1000000
        f1 = open(fastq_file, 'rU')
        #    l.acquire()
        out_file = str(i)
        for line in f1:
            self.count_read += 1
            header_1 = line.rstrip()
            if not header_1.startswith('@'):
                print(colored("Error: Sequences are not in fastq format\n", "red"))
                quit()
            read_seq = next(f1).rstrip()
            self.read_seq_orig_len_1 = len(read_seq)
            header_2 = next(f1).rstrip()
            read_qual = next(f1).rstrip()
            #   Count before filtering
            self.tot_A_1 += read_seq.count('A')
            self.tot_T_1 += read_seq.count('T')
            self.tot_G_1 += read_seq.count('G')
            self.tot_C_1 += read_seq.count('C')
            self.tot_N_1 += read_seq.count('N')
            if read_seq.count('N') > 0:
                # Reads containing N bases before filtering
                self.n_cont_read_1 += 1
            max_read_len_1, min_read_len_1 = self.find_len(read_seq, max_read_len_1, min_read_len_1)
            self.tot_len_1 += self.read_seq_orig_len_1
            GC = read_seq.count('G')+read_seq.count('C')
            self.tot_gc_1 += GC
            if self.qual_format == 1 or self.qual_format == 3:
                qual_fact = 33
                self.process_qual(read_qual, GC, qual_fact)
            elif self.qual_format == 2:
                qual_fact = 64
                self.process_qual(read_qual, GC, qual_fact)
            if len(read_seq) > int(self.min_size):
                if (float(read_seq.count('N')*100)/len(read_seq)) < int(self.n_base):
                    if self.adapter:
                        read_seq, read_qual = self.adapter_trim(read_seq, read_qual, self.Per)
                        if self.Trim:
                            read_seq, read_qual = self.qual_trim(read_seq, read_qual, self.win_size)
                            self.output_filter_data(header_1, header_2, read_seq, read_qual, out_file)
                        else:
                            read_seq, read_qual = self.qual_filter(read_seq, read_qual)
                            self.output_filter_data(header_1, header_2, read_seq, read_qual, out_file)
                    else:
                        if self.Trim:
                            read_seq, read_qual = self.qual_trim(read_seq, read_qual, self.win_size)
                            self.output_filter_data(header_1, header_2, read_seq, read_qual, out_file)
                        else:
                            read_seq, read_qual = self.qual_filter(read_seq, read_qual)
                            self.output_filter_data(header_1, header_2, read_seq, read_qual, out_file)
                else:
                    #   Reads containing more than n_base % (given)
                    self.n_read_ct_1 += 1
            else:
                self.short_read_ct_1 += 1
#       out_file.close()
        f1.close()
#       l.release()

        stat_file = open('StatTemp.txt', 'a')
        self.hash_gc_1 = collections.OrderedDict(sorted(self.hash_gc_1.items()))
        self.hash_gc_1_list = self.hash_gc_1.values()
        self.hash_gc_1a = collections.OrderedDict(sorted(self.hash_gc_1a.items()))
        self.hash_gc_1a_list = self.hash_gc_1a.values()
        self.hash_qual_1 = collections.OrderedDict(sorted(self.hash_qual_1.items()))
        self.hash_qual_1_list = self.hash_qual_1.values()
        self.hash_qual_1a = collections.OrderedDict(sorted(self.hash_qual_1a.items()))
        self.hash_qual_1a_list = self.hash_qual_1a.values()

        stat_file.write(str(self.count_read)+'\t'+str(self.trim_read_1)+'\t'+str(self.short_read_ct_1)+'\t'+'NA\t'
                       +str(self.tot_gc_1)+'\t'+str(self.tot_len_1)+'\t')    #6
        stat_file.write('NA\t'+str(self.tot_gc_a1)+'\t'+'NA\t'+str(self.tot_len_a1)+'\t'+'NA\t'+str(self.count_read_a)+
                       '\t'+str(self.n_read_ct_1)+'\t'+'NA\t'+str(self.ad_trim_1)+'\t')  #15
        stat_file.write('NA\t'+'NA\t'+str(self.fail_qual_1)+'\t'+'NA\t'+str(self.tot_A_1)+'\t'+'NA\t'+str(self.tot_T_1)+'\t'+
                       'NA\t'+str(self.tot_G_1)+'\t'+'NA\t')   #25
        stat_file.write(str(self.tot_C_1)+'\t'+'NA\t'+str(self.tot_N_1)+'\t'+'NA\t'+str(self.tot_Aa_1)+'\t'+'NA\t'+
                       str(self.tot_Ta_1)+'\t'+'NA\t'+str(self.tot_Ga_1)+'\t'+'NA\t')   #35
        stat_file.write(str(self.tot_Ca_1)+'\t'+'NA\t'+str(self.tot_Na_1)+'\t'+'NA\t'+'NA'+'\t'+'NA'+'\t'+'NA'+'\t'+
                       'NA'+'\t') #43
        stat_file.write(str(max_read_len_1)+'\t'+str(min_read_len_1)+'\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+
                       str(self.max_trim_len_1)+'\t') #52
        stat_file.write(str(self.min_trim_len_1)+'\t'+'NA\t'+str(self.tot_trim_len_1)+'\t'+'NA\t'+
                       '\t'.join(str(v) for v in self.hash_gc_1_list)+'\t')  #66
        stat_file.write('NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t')    #77
        stat_file.write('\t'.join(str(v) for v in self.hash_gc_1a_list)+'\t')   #87
        stat_file.write('NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t'+'NA\t')    #98
        stat_file.write(str(self.tot_qual_1_sum)+'\t'+str(self.tot_qual_1a_sum)+'\t'+'NA\t'+'NA\t'+str(self.n_cont_read_1)+
                       '\t'+str(self.n_cont_read_1a)+'\t'+str(self.trim_qual_ct_1)+'\t') #105
        stat_file.write('\t'.join(str(v) for v in self.hash_qual_1_list)+'\t')  #126
        stat_file.write('\t'.join(str(v) for v in self.hash_qual_1a_list)+'\n') #147
        stat_file.close()

    def find_len(self, seq, max_read_len_1L, min_read_len_1L):
        if len(seq) >= max_read_len_1L:
            max_read_len_1L = len(seq)
        if len(seq) < min_read_len_1L:
            min_read_len_1L = len(seq)
        return max_read_len_1L, min_read_len_1L

    def process_qual(self, read_qual_l, GCL, qual_fact_l):
        read_qual_l = list(read_qual_l)
        read_qual_num = map(ord, read_qual_l)
        #   For illumina 1.3+
        read_qual_num = [ele-qual_fact_l for ele in read_qual_num]
        self.tot_qual_1_sum += sum(read_qual_num)
        avg_qual = np.mean(read_qual_num)
        gc_cont = float((GCL*100)/len(read_qual_l))
        Qual1 = [ele <= 10 for ele in read_qual_num]
        Qual2 = [20 <= ele > 10 for ele in read_qual_num]
        Qual3 = [30 <= ele > 20 for ele in read_qual_num]
        Qual4 = [42 <= ele > 30 for ele in read_qual_num]
        #   self.TotQual1 += len(Qual1)
        #   self.TotQual2 += len(Qual2)
        #   self.TotQual3 += len(Qual3)
        #   self.TotQual4 += len(Qual4)
        #   histogram bins
        #   print read_qual_l
        if avg_qual == 0:
            avg_qual = 0.1
        self.hash_qual_1[int(math.ceil(avg_qual/2))*2] += 1
        if gc_cont == 0:
            gc_cont = 0.1
        self.hash_gc_1[int(math.ceil(gc_cont/10))*10] += 1

    def adapter_trim(self, read_seq_l, read_qual_l, PerL):
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
                if (float(match)/float(len(ad))) >= float(PerL):
                    del read_seq_l1[iter:iter+len_ad]
                    del read_qual_l[iter:iter+len_ad]
                    del read_seq_l[iter:iter+len_ad]
                    iter = -1
                    len2 = len(read_seq_l)
#                    read_seq_l = read_seq_l1
                iter += 1
                match = 0
            read_seq_l = ''.join(read_seq_l1)
            read_qual_l = ''.join(read_qual_l)
        if len(read_seq_l) < self.read_seq_orig_len_1:
            self.ad_trim_1 += 1
        return read_seq_l, read_qual_l

    def qual_trim(self, read_seq_l, read_qual_l, win_sizeL):
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

        while iter < (len(read_qual_num_l)-win_sizeL)+1:
            trim = read_qual_num_l[:win_sizeL]
            if np.mean(trim) <= self.qual_thresh:
                read_qual_num_l = read_qual_num_l[win_sizeL:]
                read_seq_l = read_seq_l[win_sizeL:]
                read_qual_l = read_qual_l[win_sizeL:]
                iter -= 1
            else:
                temp_read_qual_num_l.extend(read_qual_num_l[:win_sizeL])
                read_qual_num_l = read_qual_num_l[win_sizeL:]
                temp_read_seq_l.extend(read_seq_l[:win_sizeL])
                read_seq_l = read_seq_l[win_sizeL:]
                temp_read_qual_l.extend(read_qual_l[:win_sizeL])
                read_qual_l = read_qual_l[win_sizeL:]
                iter -= 1
            iter += 1
        temp_read_qual_num_l.extend(read_qual_num_l)
        temp_read_seq_l.extend(read_seq_l)
        temp_read_qual_l.extend(read_qual_l)
        if len(temp_read_seq_l) > 0 and np.mean(temp_read_qual_num_l) > self.qual_thresh:
            if len(temp_read_seq_l) < lenb:
                self.trim_qual_ct_1 += 1
            return ''.join(temp_read_seq_l), ''.join(temp_read_qual_l)
        else:
            return None, None

    def output_filter_data(self, header_1_l, header_2_l, read_seq_l, read_qual_l, out_file_l):
        self.raw_out_dir1 = self.pathname+'/'+'raw_out_1'
        os.chdir(self.raw_out_dir1)
        out_file_l = open(out_file_l, 'a')
        if out_file_l and read_seq_l and len(read_seq_l) >= int(self.min_len) and self.out_fmt == 'fastq':
            self.count_read_a += 1
            out_file_l.write(header_1_l+'\n'+read_seq_l+'\n'+header_2_l+'\n'+read_qual_l+'\n')
            self.output_filter_data_sub(read_seq_l, read_qual_l)
        elif out_file_l and read_seq_l and len(read_seq_l) > self.min_len and self.out_fmt == 'fasta':
            self.count_read_a += 1
            out_file_l.write('>'+header_1_l+'\n'+read_seq_l+'\n')
            self.output_filter_data_sub(read_seq_l, read_qual_l)

    def output_filter_data_sub(self, read_seq_l, read_qual_l):
        self.tot_Aa_1 += read_seq_l.count('A')
        self.tot_Ta_1 += read_seq_l.count('T')
        self.tot_Ga_1 += read_seq_l.count('G')
        self.tot_Ca_1 += read_seq_l.count('C')
        self.tot_Na_1 += read_seq_l.count('N')
        if read_seq_l.count('N') > 0:
            #   Reads containing N bases after filtering
            self.n_cont_read_1a += 1
        self.tot_len_a1 += len(read_seq_l)
        self.tot_gc_a1 += read_seq_l.count('G')+read_seq_l.count('C')
        gc_cont = float(((read_seq_l.count('G')+read_seq_l.count('C'))*100)/len(read_qual_l))
        if gc_cont == 0:
            gc_cont = 0.1
        self.hash_gc_1a[int(math.ceil(gc_cont/10))*10] += 1
        if len(read_seq_l) < self.read_seq_orig_len_1:
            self.trim_read_1 += 1
            self.tot_trim_len_1 += len(read_seq_l)
        self.max_trim_len_1, self.min_trim_len_1 = self.find_len(read_seq_l, self.max_trim_len_1, self.min_trim_len_1)
        read_qual_l = list(read_qual_l)
        read_qual_num_l = map(ord, read_qual_l)
        if self.qual_format == 1 or self.qual_format == 3:
            read_qual_num_l = [ele-33 for ele in read_qual_num_l]
        else:
            read_qual_num_l = [ele-64 for ele in read_qual_num_l]
        avg_qual1a = np.mean(read_qual_num_l)
        self.hash_qual_1a[int(math.ceil(avg_qual1a/2))*2] += 1
        self.tot_qual_1a_sum += sum(read_qual_num_l)

    def qual_filter(self, read_seq_l, read_qual_l):
        read_qual_l = list(read_qual_l)
        read_qual_num_l = map(ord, read_qual_l)
        if self.qual_format == 1 or self.qual_format == 3:
            read_qual_num_l = [ele-33 for ele in read_qual_num_l]
            if np.mean(read_qual_num_l) >= self.qual_thresh:
                return read_seq_l, ''.join(read_qual_l)
            else:
                self.fail_qual_1 += 1
                return None, None
        elif self.qual_format == 2:
            read_qual_num_l = [ele-64 for ele in read_qual_num_l]
            if np.mean(read_qual_num_l) >= self.qual_thresh:
                return read_seq_l, ''.join(read_qual_l)
            else:
                self.fail_qual_1 += 1
                return None, None

    def detect_fastq_variant(self):
        Count = 0
        Check = []
        #File = open(self.pathname+'/'+self.file_1, 'rU')
        File = open(self.file_1, 'rU')

        for line in File:
            id =line.rstrip()
            if not id.startswith('@'):
                print(colored("Error: Sequences are not in fastq format\n", "red"))
                sys.exit()
            seq = next(File).rstrip()
            #   input parameter check
            if len(seq) <= int(self.min_size):
                print(colored("Error: The minimum length is greater or equal to actual sequence length", "red"))
                sys.exit()
            next(File)
            Asc = next(File).rstrip()
            asc_list = list(Asc)
            asc_list = list(map(ord, asc_list))
            Min = min(asc_list)
            Max = max(asc_list)
            Check.append(Min)
            Check.append(Max)
            Count += 1
            if Count == 40000:
                break
        File.close()
        Min = min(Check)
        Max = max(Check)
        if 64 > Min >= 33 and Max == 74:
            return 1
        elif Min >= 64 and 74 < Max <= 104:
            return 2
        elif 64 > Min >= 33 and Max <= 73:
            return 3

    def get_hash(self, start, end):
        d = {}
        for v in range(start, end, start):
            d[v] = 0
        return d

    def file_parts(self, iterable, n):
        iterable = iter(iterable)
        while True:
            yield chain([next(iterable)], islice(iterable, n-1))

    @staticmethod
    def usage():
        print("\nusage: srap -filter-s [options] -a fastq_file")

        print("\nOptions:")
        print("  "+"-a, --p1 INT {:>34}".format("Input file (.fastq, .fq)"))
        print("  "+"-d, --msz INT {:>60}".format("filter the reads which are lesser than minimum size"))
        print("  "+"-c, --qfmt INT {:>29}".format("Quality value format "))
        print("{:>40}".format("1= Illumina 1.8"))
        print("{:>40}".format("2= Illumina 1.3"))
        print("{:>35}".format("3= Sanger\n"))
        print("  "+"-e, --nb INT {:>51}".format("filter the reads containing given % of N "))
        print("  "+"-f, --adp STRING {:>81}".format("Trim the adapter sequence and truncate the read "
                                                    "sequence [adapter sequence]"))
        print("  "+"-g, --per FLOAT {:>124}".format("Truncate the read sequence if it matches to adapter sequence "
                                                    "equal or more than given percent (0.0-1.0) [default=0.9]"))
        print("  "+"-i, --qthr INT {:>113}".format("Filter the read sequence if average quality of "
                                                   "bases in reads is lower than threshold (1-40) [default:20]"))
#        print "  "+"-j, --mqual INT {:>122}".format("filter the reads if given percentage of bases having lower "
#                                                    "quality values lower than threshold (1-100) [default:20]")
        print("  "+"-n, --trim BOOLEAN {:>143}".format("If trim option set to true, the reads with low "
                                                       "quality (as defined by option --qthr) will be trimmed instead"
                                                       " of discarding [default: False]"))
        print("  "+"-p, --wsz INT {:>118}".format("The window size for trimming (5\'->3\') the reads. This option"
                                                  " should always set when -trim option is defined. "))
        print("{:>101}".format("The recommended windowsize for best result should be between 2-5 [default:5]"))
        print("  "+"-r, --mlk INT {:>65}".format("The minimum length of the reads to retain after trimming"))
        print("  "+"-q, --cpu {:>38}".format("Number of CPU [default:2]"))
        print("  "+"-m, --ofmt {:>60}".format("Output file format (fastq/fasta) [default:fastq]"))
        print("  "+"-v, --no-vis BOOLEAN {:>51}".format("No figures will be produced [yes|no] [default:no]"))
        print("  "+"-h, --help {:>37}".format("Print this help message\n\n"))


if __name__ == '__main__':
    ob = FilterSingle()
    try:
        ob.filter_single()
    except KeyboardInterrupt:
        print(colored("\nProgram is terminated by user", "red"))
        sys.exit(1)

