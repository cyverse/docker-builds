#!/usr/bin/env python

###############################################################################
#                                                                             #
#    BatchBowtie                                                              #
#                                                                             #
#    A wrapper script for Bowtie2, that runs Bowtie and Samtools on a         #
#    directory of read files.                                                 #
#                                                                             #
#    Copyright (C) Benjamin Bolduc                                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This library is free software; you can redistribute it and/or            #
#    modify it under the terms of the GNU Lesser General Public               #
#    License as published by the Free Software Foundation; either             #
#    version 3.0 of the License, or (at your option) any later version.       #
#                                                                             #
#    This library is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        #
#    Lesser General Public License for more details.                          #
#                                                                             #
#    You should have received a copy of the GNU Lesser General Public         #
#    License along with this library.                                         #
#                                                                             #
###############################################################################

__author__ = "Ben Bolduc"
__copyright__ = "Copyright 2016"
__credits__ = ["Ben Bolduc"]
__license__ = "LGPLv3"
__maintainer__ = "Ben Bolduc"
__email__ = "bolduc.10@osu.edu"
__status__ = "Development"

import os
import sys
import subprocess
import argparse
import itertools
from pprint import pprint
from natsort import natsorted
from Bio import SeqIO
import Levenshtein

parser = argparse.ArgumentParser(description="The script essentially wraps bowtie2 for aligning an entire directory of "
                                             "paired or unpaired reads against the same reference sequence collection."
                                             " This new version now handles paired and unpaired reads simultaneously. "
                                             "Interleaved files (files including both forward and reverse reads) are "
                                             "separated for Bowtie.",
                                 formatter_class=argparse.RawTextHelpFormatter)

inputs = parser.add_argument_group('Required Inputs and Parameters')

inputs.add_argument('--reads', '--reads-dir', dest='reads_dir', metavar='DIRECTORY', default='/workDir',
                    help="Directory containing the reads files to be aligned. Selecting a single file instead of a "
                         "directory will result in only one file being read.")

inputs.add_argument('-d', '--db', '--input-db', dest='input_db', metavar='FILENAME',
                    help="The NAME of a FASTA-formatted sequence file containing sequences/references reads to be "
                         "aligned against. For example, if the filename is 'contigs.fasta' then '--db contigs.fasta'")

inputs.add_argument('-f', '--fmt', '--input-format', dest='input_fmt', choices=['fastq', 'fasta', 'fq'],
                    help="File format of reads to be aligned. Compressed files (*.gz, *.tar.gz) will be automatically"
                    "recognized.")

inputs.add_argument('-s', '--interleaved', dest='separate', action='store_true',
                    help="If enabled, will treat each file as INTERLEAVED and process them as paired files. "
                         "Interleaved reads are assumed to be in the format F,R,F,R.")

genOpts = parser.add_argument_group('General Options')  # Bowtie2 and Samtools no longer specified

genOpts.add_argument('--read-types', dest='read_types', choices=['paired', 'unpaired', 'mixed'], default='mixed',
                     help="Whether or not reads are paired, unpaired or mixed.")

genOpts.add_argument('--dist', dest='distance', type=int, default=1,
                     help="Levenshtein distance between filenames, used to determine pairing/grouping of paired and"
                          " unpaired files. So OSD101_R1 and OSD101_R2 have a distance of 1, OSD101_R1_paired and "
                          "OSD101_R2_unpaired have a distance of 2. For many cases (like using all paired reads), the "
                          "default is fine. It will correctly identify all the read 'pairings.' For mixing paired and "
                          "unpaired files, it is suggested to use 2-3. That way the script can identify AND GROUP all "
                          "the paired and unpaired reads from a single experiment into a single group. Basically, "
                          "naming within a group should be closer than between groups. Naturally sorting files ASSUMES "
                          "that unpaired read files will ALWAYS follow their paired version.")

genOpts.add_argument('-x', '--exclude', dest='filter', action='append', default=['.py'],
                     help="A list of strings (with -x for each item) to filter out files within the input directory. "
                          "This is useful when there are paired, unpaired, QCd and non-QCd versions of read files "
                          "existing in the soure directory - and the user only wants finalized, paired reads. For "
                          "example, one can use '-x unpaired -x singles' to remove all files containing the words"
                          " 'unpaired' and 'singles.' The string matching is not case-sensitive.")

genOpts.add_argument('-u', '--unpair-terms', dest='unpair_term', action='append', default=['unpair'],
                     help="A list of strings (with -u for each item) to help identify unpaired read files (if in the "
                          "input). By default, unpaired files are identified by their natural sort order against "
                          "'predicted' paired-end files. For example, if all unpaired read files have 'unpaired' in "
                          "their name, you can use '-u unpaired' The string matching is not case-sensitive.")

genOpts.add_argument('-p', '--pair-terms', dest='unpair_term', action='append', default=['pair'],
                     help="A list of strings (with -p for each item) to help identify paired read files (if in the "
                          "input). By default, paired files are identified by their natural sort order against "
                          "'predicted' unpaired-end files. For example, if all paired read files have 'paired' in "
                          "their name, you can use '-p paired' This is applied AFTER unpaired reads are filtered out."
                          " The string matching is not case-sensitive.")

genOpts.add_argument('-k', '--keep-sam', dest='keep_sam', action='store_true',
                     help="If enabled, SAM files will be preserved during BAM file generation. Without this option, "
                          "THERE WILL BE NO SAM FILES.")

genOpts.add_argument('-m', '--merge-results', dest='merge_output', action='store_true',
                     help="If enabled, and if multiple reads are supplied, combine the bowtie2 results into one file.")

genOpts.add_argument('--merge-name', dest='merge_name', metavar='FILENAME', default='bowtie2-run.sam',
                     help="Filename to use for merged results. This WILL NOT be used if --merge-results isn't used.")

genOpts.add_argument('-z', '--remove-tmp', dest='remove_tmp', action='store_true',
                     help="Remove read files used to generate BAM files. This includes the intermediary bowtie2"
                          " database files and the fastq (or other format) read files. This option should only be "
                          "enabled if the system COPIES source files into another location, as it standard within the "
                          "Cyverse cyberinfrastructure.")

genOpts.add_argument('-l', '--log-file', dest='log_fn', metavar='FILENAME', default='bowtie2-read-mapping.log',
                    help="Log file name")

bowtie2Opts = parser.add_argument_group('Bowtie2 Alignment Options')

bowtie2Opts.add_argument('--alignment_type', dest='align_type', choices=['end-to-end', 'local'], default='end-to-end',
                         help="Whether the entire read must align (end-to-end) or only a local region (local).")

bowtie2Opts.add_argument('--end-to-end_presets', dest='global_presets', metavar='STRING',
                         choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
                         default="sensitive", help="Presets for end-to-end alignments:\n"
                         "very-fast: -D 5 -R 1 -N 0 -L 22 -i S,0,2.50\n"
                         "fast: -D 10 -R 2 -N 0 -L 22 -i S,0,2.50\n"
                         "sensitive: -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)\n"
                         "very-sensitive: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50")

bowtie2Opts.add_argument('--local_presets', dest='local_presets', metavar='STRING',
                         choices=['very-fast-local', 'fast-local', 'sensitive-local', 'very-sensitive-local'],
                         default='sensitive-local',
                         help="Presets for local alignments:\n"
                         "very-fast-local: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00\n"
                         "fast-local: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75\n"
                         "sensitive-local: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)\n"
                         "very-sensitive-local: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50")

bowtie2Opts.add_argument('--non-deterministic', dest='non_deterministic', action='store_true',
                     help="Bowtie 2 will use the current time to re-initialize the pseudo-random number generator. "
                          "Useful when the input consists of many identical reads.")

bowtie2Opts.add_argument('--trim5', dest='trim5', metavar='INT', type=int, default=0,
                         help="Trim X bases from 5'/left end of reads.")
bowtie2Opts.add_argument('--trim3', dest='trim3', metavar='INT', type=int, default=0,
                         help="Trim X bases from 3'/right end of reads.")
bowtie2Opts.add_argument('-I', '--minins', dest='minins', metavar='INT', type=int, default=0,
                         help="Minimum fragment length.")
bowtie2Opts.add_argument('-X', '--maxins', dest='maxins', metavar='INT', type=int, default=2000,
                         help="Maximum fragment length.")

results = parser.parse_args()


def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)

# Ensure there's a DB file
if not results.input_db:
    error("An input bowtie2 database file or fasta file is required.")

preset = False
if results.align_type == "end-to-end":
    preset = results.global_presets
if results.align_type == "local":
    preset = results.local_presets

if results.input_fmt not in ['fasta', 'fastq', 'fq']:
    error('ERROR: Input format must be either fasta or fastq formatted')

# TODO Write additional code to handle interleaved + non-paired data
if (results.read_types == 'mixed') and results.separate:
    error('ERROR: cannot use mixed read types and interleaved. NOT IMPLEMENTED YET.')

if (results.read_types == 'unpaired') and results.separate:
    error('ERROR: cannot use unpaired read types and interleaved. NOT IMPLEMENTED YET')


def execute(command, logfile):

    logfile.write('Executing {}'.format(command) + os.linesep)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()

    logfile.write(stderr + os.linesep)


def file_finder(rootdir, searchstring):
    found_list = []

    for root, dirs, files in os.walk(rootdir):
        # found_list = [os.path.join(root, x) for x in files if searchstring in x]
        for x in files:
            if searchstring in x:
                found_list.append(os.path.join(root, x))

    return found_list


def group(lst, n):
    """http://code.activestate.com/recipes/303060-group-a-list-into-sequential-n-tuples/"""
    # Why aren't I using grouper from itertools recipe?
    return list(zip(*[lst[i::n] for i in range(n)]))


def split_interleaved(reads_lst, logfile):

    read_type = results.readType
    if results.readType in ['fastq', 'fq']:
        read_type = 'fastq'
    # Doesn't matter if it's fasta, because read_type will inherit that

    paired_reads = []
    for reads in reads_lst:
        logfile.write('Splitting {}'.format(reads) + os.linesep)
        reads1_fn = reads.rsplit('.')[0] + '_R1.{}'.format(read_type)
        reads2_fn = reads.rsplit('.')[0] + '_R2.{}'.format(read_type)
        reads1 = []
        reads2 = []

        with open(reads, 'rU') as reads_fh:
            for read1, read2 in group(SeqIO.parse(reads_fh, read_type), 2):
                reads1.append(read1)
                reads2.append(read2)

        with open(reads1_fn, 'w') as reads1_fh:
            SeqIO.write(reads1, reads1_fh, read_type)
        with open(reads2_fn, 'w') as reads2_fh:
            SeqIO.write(reads2, reads2_fh, read_type)

        logfile.write('Split files:\n{}\n'.format(reads1_fn, reads2_fn) + os.linesep)
        paired_reads.append((reads1_fn, reads2_fn))

    return paired_reads


def group_reads(reads_lst, logfile):

    natsorted_list = natsorted(reads_lst)

    reads_groups = {}

    for i in range(len(natsorted_list) - 1):
        string1 = natsorted_list[i]
        string2 = natsorted_list[i + 1]
        dist = Levenshtein.distance(string1, string2)

        if dist <= results.distance:

            if not reads_groups:  # Initialize
                reads_groups[0] = {string1, string2}

            updated = False
            for group_num, group_set in reads_groups.items():

                if not group_set.isdisjoint({string1, string2}):  # If there is something in common, it's NOT empty
                    reads_groups[group_num].update({string1, string2})
                    updated = True

            if not updated:  # If NONE of the elements were found in ANY of the sets
                reads_groups[len(reads_groups.keys())] = {string1, string2}

    organized_reads = {}

    # Identify which are paired and unpaired...
    for group_num, group_set in reads_groups.items():
        # Set generation occurred during incrementing through naturally sorted list... so there shouldn't be any change
        re_natural_sort = natsorted(group_set)

        if group_num not in organized_reads:
            organized_reads[group_num] = {}

        if results.read_types == 'mixed':

            # If there's 4 in a "set", then it *should* be paired, unpaired, paired, unpaired
            if len(re_natural_sort) == 4:

                if 'paired' not in organized_reads[group_num]:
                    organized_reads[group_num]['paired'] = re_natural_sort[::2]

                if 'unpaired' not in organized_reads[group_num]:
                    organized_reads[group_num]['unpaired'] = re_natural_sort[1::2]

            if len(re_natural_sort) == 2:
                logfile.write(
                    'WARNING: "Mixed" specified, yet there are only 2 files identified. Processing as if the read files'
                    'are PAIRED. If they are both unpaired, specify "unpaired" as the read-type. If the grouping is not'
                    'correct, try adjusting the distance.' + os.linesep)

                if 'paired' not in organized_reads[group_num]:
                    organized_reads[group_num]['paired'] = re_natural_sort

            if len(re_natural_sort) == 3:  # Wonderful...

                unpaired = []
                # See if user provided any useful search terms
                for read_fn in re_natural_sort:
                    for filterStr in results.unpair_term:
                        if filterStr in read_fn:

                            if 'unpaired' not in organized_reads[group_num]:
                                organized_reads[group_num]['unpaired'].append(read_fn)

                            re_natural_sort.remove(read_fn)

                            if len(re_natural_sort) == 2:  # Now that the oddball was removed

                                if 'paired' not in organized_reads[group_num]:
                                    organized_reads[group_num]['paired'] = re_natural_sort

                                else:
                                    logfile.write(
                                        'ERROR: "Mixed" specified yet 1 or 3 files remaining that cannot be classified '
                                        'as paired or unpaired reads. If the grouping is not correct, try adjusting the'
                                        ' distance OR adding unpaired terms to help assist the program in identifying '
                                        'unpaired/paired reads.')
                                    sys.exit(1)

            if len(re_natural_sort) == 1:
                logfile.write(
                    'WARNING: "Mixed" specified, yet there is only 1 file identified. Processing as if the read file is'
                    ' UNPAIRED. If this grouping is not correct, try adjusting the distance.')

                if 'unpaired' not in organized_reads[group_num]:
                    organized_reads[group_num]['unpaired'] = re_natural_sort

        if results.read_types == 'paired':
            # Ensure pairs are only in PAIRS
            if len(re_natural_sort) != 2:
                if not results.separate:
                    logfile.write(
                        'ERROR: "Paired" specified, yet read pairs arent in groups of 2. Perhaps change the distance?')
                    sys.exit(1)
                else:
                    logfile.write(
                        'NOTICE: "Paired" specified, Identified interleaved files.')
                    if 'paired' not in organized_reads[group_num]:
                        organized_reads[group_num]['paired'] = re_natural_sort

            if len(re_natural_sort) == 2:

                if 'paired' not in organized_reads[group_num]:
                    organized_reads[group_num]['paired'] = re_natural_sort

        if results.read_types == 'unpaired':

            if 'unpaired' not in organized_reads[group_num]:
                organized_reads[group_num]['unpaired'] = re_natural_sort

    return organized_reads


def prepare_reads(reads_dict, reads_dir, logfile):
    # Will have full paths for each of the files
    # num : {paired: [], unpaired: []}
    new_locations = reads_dict

    for group_num, group_dict in reads_dict.items():

        if 'paired' in group_dict:  # Don't bother checking if they're 0

            if len(group_dict['paired']) == 2:

                paired_list = group_dict['paired']

                paired1, paired2 = paired_list

                if all('.gz' in pair for pair in paired_list):  # If gz is in all the pairs in paired list

                    if not all('.tar.gz' in pair for pair in paired_list):  # If it's got gz AND .tar.gz

                        # Will zcat next, so need the 'new' filename - Don't want entire path - it has source dir!
                        input1_fn = os.path.join(reads_dir, os.path.basename(paired1).rsplit('.', 1)[0])
                        input2_fn = os.path.join(reads_dir, os.path.basename(paired2).rsplit('.', 1)[0])

                        if os.path.exists(input1_fn) and os.path.exists(input2_fn):
                            logfile.write('File already exists at designated location' + os.linesep)

                            new_locations[group_num]['paired'] = [input1_fn, input2_fn]
                            continue

                        cmd = ''
                        if os.name == 'posix':  # Oh the time spent figuring this out
                            cmd = 'zcat < {0} > {1} && zcat < {2} > {3}'.format(paired1, input1_fn, paired2, input2_fn)
                        else:
                            cmd = 'zcat {0} > {1} && zcat {2} > {3}'.format(paired1, input1_fn, paired2, input2_fn)

                        execute(cmd, logfile)

                        new_locations[group_num]['paired'] = [input1_fn, input2_fn]

                    else:  # Then it only has .gz

                        if all('.tar.gz' in pair for pair in paired_list):
                            # Will zcat next, so need the 'new' filename - Don't want entire path - it has source dir!
                            input1_fn = os.path.join(reads_dir, os.path.basename(paired1).rsplit('.', 2)[0])
                            input2_fn = os.path.join(reads_dir, os.path.basename(paired2).rsplit('.', 2)[0])

                            if os.path.exists(input1_fn) and (os.path.exists(input2_fn)):
                                logfile.write('File already exists at designated location' + os.linesep)
                                new_locations[group_num]['paired'] = [input1_fn, input2_fn]
                                continue

                            cmd = 'tar xf {0} -C {1} && xf {2} -C {1}'.format(paired1, reads_dir, paired2)
                            execute(cmd, logfile)

                            new_locations[group_num]['paired'] = [input1_fn, input2_fn]

                        else:
                            logfile.write("ERROR: There's .gz in the filenames and 'not' .tar.gz, but not all the "
                                          ".tar.gz files are... in all of them?")
                            sys.exit(1)

                else:  # If not .gz files then assuming it's uncompressed... so this will fail for .zip, .bz2, etc...

                    input1_fn = os.path.join(reads_dir, os.path.basename(paired1))
                    input2_fn = os.path.join(reads_dir, os.path.basename(paired2))

                    if os.path.exists(input1_fn) and (os.path.exists(input2_fn)):
                        logfile.write('File already exists at designated location' + os.linesep)
                        new_locations[group_num]['paired'] = [input1_fn, input2_fn]
                        continue

                    cmd = 'cp {0} {1} && cp {2} {1}'.format(paired1, reads_dir, paired2)
                    execute(cmd, logfile)

                    new_locations[group_num]['paired'] = [input1_fn, input2_fn]

            elif (len(group_dict['paired']) == 1) and results.separate:
                logfile.write('Interleaved files enabled....' + os.linesep)
                output_fn = ''
                interleaved_fn = group_dict[group_num]['paired'][0]  # It's a list!

                if '.tar.gz' in group_dict[group_num]['paired']:
                    output_fn = os.path.join(reads_dir, os.path.basename(interleaved_fn).rsplit('.', 2)[0])
                    cmd = 'tar xf {0} -C {1}'.format(interleaved_fn, reads_dir)
                    execute(cmd, logfile)

                # Once decompression is finished, can run through biopython
                paired_reads = split_interleaved([output_fn],
                                                 logfile)  # Returns list of tuples... well, a list containing 1 tuple
                print('There\'s paired_reads for splitting interleaved?')
                print(paired_reads[0], paired_reads[1])
                new_locations[group_num]['paired'] = [paired_reads[0], paired_reads[1]]

            else:
                logfile.write('ERROR: During preparing reads for bowtie. Length of paired reads was greater than 2.')
                sys.exit(1)

        if 'unpaired' in group_dict:

            unpaired_list = group_dict['unpaired']

            for pos, unpaired in enumerate(unpaired_list):  # Unlike dict, don't know 'location' of read in list

                cmd = ''

                if '.tar.gz' in unpaired:

                    unpaired_fn = os.path.join(reads_dir, os.path.basename(unpaired.rsplit('.', 2)[0]))  # Output file

                    if os.path.exists(unpaired_fn):
                        logfile.write('{} already exists at designated location'.format(unpaired_fn) + os.linesep)
                        new_locations[group_num]['unpaired'][pos] = unpaired_fn
                        continue

                    cmd = 'tar xf {0} -C {1}'.format(unpaired, reads_dir)
                    execute(cmd, logfile)

                    new_locations[group_num]['unpaired'][pos] = unpaired_fn

                elif '.gz' in unpaired:

                    unpaired_fn = os.path.join(reads_dir, os.path.basename(unpaired.rsplit('.', 1)[0]))  # Output file

                    if os.path.exists(unpaired_fn):
                        logfile.write('{} already exists at designated location'.format(unpaired_fn) + os.linesep)
                        new_locations[group_num]['unpaired'][pos] = unpaired_fn
                        continue

                    cmd = ''

                    if os.name == 'posix':  # Oh the time spent figuring this out
                        cmd = 'zcat < {0} > {1}'.format(unpaired, unpaired_fn)
                    else:
                        cmd = 'zcat {0} > {1}'.format(unpaired, unpaired_fn)

                    execute(cmd, logfile)
                    new_locations[group_num]['unpaired'][pos] = unpaired_fn

                else:  # Treat as any other reads file

                    unpaired_fn = os.path.join(reads_dir, os.path.basename(unpaired))

                    if os.path.exists(unpaired_fn):
                        logfile.write('{} already exists at designated location'.format(unpaired_fn) + os.linesep)
                        new_locations[group_num]['unpaired'][pos] = unpaired_fn
                        continue

                    cmd = 'cp {0} {1}}'.format(unpaired, reads_dir)
                    execute(cmd, logfile)

                    new_locations[group_num]['unpaired'][pos] = unpaired_fn

    return new_locations


def prepare_bowtie_db(fasta, db_dir, logfile):

    db_list = os.listdir(db_dir)  # No full paths

    if '.DS_Store' in db_list:
        db_list.remove('.DS_Store')

    bt2_db_fasta = os.path.join(db_dir, os.path.basename(fasta))
    bt2_db_base = os.path.join(db_dir, bt2_db_fasta.rsplit('.', 1)[0])

    if os.path.basename(fasta) in db_list:  # Fasta file already copied over
        logfile.write('WARNING: DB fasta file already present in bowtie2 database directory. Skipping' + os.linesep)
        db_list.remove(os.path.basename(fasta))

    else:
        # Move fasta file to bowtie2-db and create DB
        copy_fasta_cmd = 'cp {} {}'.format(fasta, db_dir)
        execute(copy_fasta_cmd, logfile)

    if db_list:

        if all('.bt2' in fn for fn in db_list):  # Understand WHY all(empty list) = True... but really?
            logfile.write('WARNING: Bowtie2 DB files already present in bowtie2 database directory. Skipping' + os.linesep)

        else:
            bowtie_db_cmd = 'bowtie2-build -f {} {}'.format(bt2_db_fasta, bt2_db_base)
            execute(bowtie_db_cmd, logfile)

    else:
        bowtie_db_cmd = 'bowtie2-build -f {} {}'.format(bt2_db_fasta, bt2_db_base)
        execute(bowtie_db_cmd, logfile)

    return bt2_db_base


def merged_bowtie(reads_dict, bowtie2_db):
    """At this point should have a list of paired reads that will need to be aligned by bowtie2."""
    # num : {paired: [], unpaired: []}

    pairs1_list = [group_type['paired'][::2] for group_num, group_type in reads_dict.items() if 'paired' in group_type]
    pairs2_list = [group_type['paired'][1::2] for group_num, group_type in reads_dict.items() if 'paired' in group_type]

    # The list is really a list of lists (one for each group), need to flatten
    bowtie_pairs1 = ','.join(list(itertools.chain.from_iterable(pairs1_list)))
    bowtie_pairs2 = ','.join(list(itertools.chain.from_iterable(pairs2_list)))

    unpair_list = [group_type['unpaired'] for group_num, group_type in reads_dict.items() if 'unpaired' in group_type]
    bowtie_unpairs = ','.join(list(itertools.chain.from_iterable(unpair_list)))

    inFmt = False
    if 'fasta' == results.input_fmt:
        inFmt = '-f'
    if 'fastq' == results.input_fmt:
        inFmt = '-q'
    if 'fq' == results.input_fmt:
        inFmt = '-q'

    if not inFmt:
        error('No format selected during bowtie2 search.')

    input_cmd = ''

    if (len(pairs1_list) > 0) and (len(pairs2_list) > 0):
        input_cmd += ' -1 {} -2 {}'.format(bowtie_pairs1, bowtie_pairs2)
    if len(bowtie_unpairs) > 0:
        input_cmd += ' -U {}'.format(bowtie_unpairs)

    if results.non_deterministic:
        input_cmd += ' --non-deterministic {}'

    bowtie2_cmd = 'bowtie2 {} --phred33 --{} --{} -p 16 -I {} -X {} --no-unal -x {}{}'.format(
        inFmt, results.align_type, preset, results.minins, results.maxins, bowtie2_db, input_cmd)

    sam_out = results.merge_name

    return [(bowtie2_cmd, sam_out)]


def separate_bowtie(reads_dict, bowtie2_db):

    bowtie2_cmds = []

    inFmt = False
    if 'fasta' == results.input_fmt:
        inFmt = '-f'
    if 'fastq' == results.input_fmt:
        inFmt = '-q'
    if 'fq' == results.input_fmt:
        inFmt = '-q'

    if not inFmt:
        error('No format selected during bowtie2 search.')

    # num : {paired: [], unpaired: []}
    for group_num, group_type in reads_dict.items():

        input_cmd = ''
        sam_out = ''

        if 'paired' in group_type:  # Shouldn't have to double-check at this point
            input_cmd += ' -1 {} -2 {}'.format(group_type['paired'][0], group_type['paired'][1])
            sam_out = group_type['paired'][0].rsplit('.', 1)[0] + '.sam'
        if 'unpaired' in group_type:
            input_cmd += ' -U {}'.format(','.join(group_type['unpaired']))

            if not sam_out:
                sam_out = group_type['unpaired'][0].rsplit('.', 1)[0] + '.sam'

        bowtie2_cmd = 'bowtie2 {} --phred33 --{} --{} -p 16 -I {} -X {} --no-unal -x {}{}'.format(
        inFmt, results.align_type, preset, results.minins, results.maxins, bowtie2_db, input_cmd)

        if results.non_deterministic:
            bowtie2_cmd += ' --non-deterministic {}'

        bowtie2_cmds.append((bowtie2_cmd, sam_out))

    return bowtie2_cmds


def run_bowtie(cmd2run, keep_sam, logfile):

    processCall = ''

    if len(cmd2run) > 1:  # Alternatively, if sam_output

        for (bowtie2, sam_out) in cmd2run:

            bam_out = sam_out.replace('sam', 'bam')

            # Keep sam file or go directly to bam. Two ways of keeping sam file... one of which sometimes hasn't worked

            if keep_sam:
                processCall = '{0} -S {1} && samtools view -Sb -o {2} {1}'.format(bowtie2, sam_out, bam_out)

            else:
                processCall = '{0} | samtools view -Sb - > {1}'.format(bowtie2, bam_out)

            execute(processCall, logfile)

    if len(cmd2run) == 1:  # just double-checking

        bowtie2, sam_out = cmd2run[0]
        bam_out = sam_out.replace('sam', 'bam')

        if keep_sam:
            processCall = '{0} -S {1} && samtools view -Sb -o {2} {1}'.format(bowtie2, sam_out, bam_out)

        else:
            processCall = '{0} | samtools view -Sb - > {1}'.format(bowtie2, bam_out)
            #output = subprocess.check_output(processCall, shell=True)  # Trusted input only!

        execute(processCall, logfile)  # bowtie uses stderr to print out results

    if (len(cmd2run) < 1) and (len(cmd2run) != 1):  # This *should not* ever be triggered
        logfile.write('ERROR: Not >1 or ==1, but sam file specified')

if __name__ == '__main__':

    # Set up directories
    cwd = os.getcwd()
    reads_scratch = os.path.join(cwd, 'reads')
    db_scratch = os.path.join(cwd, 'bowtie2-db')

    log = open(results.log_fn, 'w', 0)  # Set buffer size to 0 to force flushing to disk

    if not os.path.exists(reads_scratch):
        os.makedirs(reads_scratch)
        log.write('Created reads directory: {}\n'.format(reads_scratch))

    if not os.path.exists(db_scratch):
        os.makedirs(db_scratch)
        log.write('Created bowtie2-db directory: {}\n'.format(db_scratch))

    log.write('Directory content for reads:' + os.linesep)
    pprint(os.listdir(reads_scratch), log)

    log.write('Directory contents for bowtie2-db directory:' + os.linesep)
    pprint(os.listdir(db_scratch), log)

    workable_files = file_finder(results.reads_dir, results.input_fmt)

    log.write('Workable Files' + os.linesep)
    pprint(workable_files, log)

    # Filter
    workable_reads = [files for files in workable_files if not any(filterStr in files for filterStr in results.filter)]
    log.write('Workable Reads' + os.linesep)
    pprint(workable_reads, log)

    # Bowtie2 DB will be in 'root' directory, whereas all the processing of the files will be in the workdir
    bt2_db_base = prepare_bowtie_db(results.input_db, db_scratch, log)
    log.write('Bowtie2 base db: {}'.format(bt2_db_base) + os.linesep)

    grouped_workable_reads = group_reads(workable_reads, log)

    log.write('Sorted files:' + os.linesep)
    pprint(grouped_workable_reads, log)

    # All reads are grouped, now just generate their file locations after decompressing, if necessary
    read_locs = prepare_reads(grouped_workable_reads, reads_scratch, log)  # Already have full paths
    log.write('Reads in their appropriate locations?' + os.linesep)
    pprint(read_locs, log)

    cmd_and_sam = []

    if results.merge_output:
        # Run merged_bowtie
        cmd_and_sam = merged_bowtie(read_locs, bt2_db_base)

    else:  # Individually run
        cmd_and_sam = separate_bowtie(read_locs, bt2_db_base)

    run_bowtie(cmd_and_sam, results.keep_sam, log)

    if results.remove_tmp:
        for files in workable_reads:
            os.remove(files)

        for files in os.listdir(db_scratch):
            os.remove(files)

        os.removedirs(db_scratch)

    log.write('Program Complete')
    log.close()
