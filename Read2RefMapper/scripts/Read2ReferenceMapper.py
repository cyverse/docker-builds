#!/usr/bin/env python

###############################################################################
#                                                                             #
#    Read2ReferenceMapper                                                     #
#                                                                             #
#    A wrapper script, written for Docker, that combines several tools to     #
#    efficiently map, parse, and filter reads against a set of reference      #
#    sequences.                                                               #
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
__credits__ = ["Ben Bolduc", "Simon Roux"]
__license__ = "LGPLv3"
__maintainer__ = "Ben Bolduc"
__email__ = "bolduc.10@osu.edu"
__status__ = "Development"

import os
import subprocess
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Don't use X-Windows backend
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from pprint import pprint
from palettable.colorbrewer.sequential import PuBuGn_9

parser = argparse.ArgumentParser(description=
                                 'A wrapper script for BamM for use on the OSC HPC. The user needs to have '
                                 'samtools and BamM in their path to run this wrapper. This also wraps Simon\'s'
                                 'filter_bam_file_coverage.py script, which filters bam files to keep only reads mapping'
                                 'to sequences coverage on at least XX % of their length')

program_group = parser.add_argument_group('Required Options')
program_group.add_argument('--dir', required=True, help='Directory with BAM files. They MUST have the BAM extension')
program_group.add_argument('--log', dest='log_fn', help='Logging file.')
program_group.add_argument('--num-threads', dest='threads', default=4, help='Threads to use for calculations.')
program_group.add_argument('--metagenome-sizes', dest='metagenome_size',
                           help='A csv-formatted file containing metagenome name and metagenome size. If provided, '
                                '(along with --coverages), will normalize coverages to the size of their metagenomes.')

cov_filter_group = parser.add_argument_group('Filter Coverage Options')
cov_filter_group.add_argument('--cov_filter', default=0, help='% length of sequences reads must map against.')

parse_group = parser.add_argument_group('BamM Parse Options')
parse_group.add_argument('--links', help='Filename to write pairing links')
parse_group.add_argument('--inserts', help='Filename to write insert distributions')
parse_group.add_argument('--coverages', default='coverage_table.csv', help='Filename to write coverage profiles')
parse_group.add_argument('--num-types', default=1, type=int, help='Number of insert/orientation types per BAM')

filter_group = parser.add_argument_group('BamM Filter Options')
filter_group.add_argument('--percent-id', dest='percent_id', type=float,
                          help='Minimum allowable percentage base identity of a mapped read')
filter_group.add_argument('--percent-aln', dest='percent_aln', type=float,
                          help='Minimum allowable percentage read bases mapped')

filter_group.add_argument('--coverage-mode', dest='coverage_mode', default='tpmean',
                          help='How to calculate coverage* (req --coverages)')
filter_group.add_argument('--cutoff-range', dest='cutoff-range',
                          help='Used to calculate upper / lower rejection cutoffs when calculating coverage')
filter_group.add_argument('--length', type=int, help='Minimum Q length')
filter_group.add_argument('--base_quality', type=int, help='Base quality threshold (Qscore)')
filter_group.add_argument('--mapping_quality', type=int, help='Mapping quality threshold')
filter_group.add_argument('--max_distance', type=int, help='Maximum allowable edit distance from query to reference')

results = parser.parse_args()


def file_finder(rootdir, searchstring):
    found_list = []

    for root, dirs, files in os.walk(rootdir):
        for x in files:
            if searchstring in x:
                found_list.append(os.path.join(root, x))
        #found_list = [os.path.join(root, x) for x in files if searchstring in x]
    return found_list


def execute(command, logfile):

    logfile.write('Executing {}'.format(command))
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()

    return stdout, stderr


def cp_files(fileList, targetDir, logfile):
    """Because it's the first stage in the process, the files must first be copied over to $TMPDIR"""
    scratch_files = []

    for files in fileList:
        copy_cmd = 'cp {} {}'.format(files, targetDir)

        cp_stdout, cp_stderr = execute(copy_cmd, logfile)
        logfile.write(cp_stderr)

        fileOut = os.path.join(targetDir, os.path.basename(files))
        scratch_files.append(fileOut)

    return scratch_files


def bamm_filter(sortedBamList, destDir, logfile):

    filteredBams = []
    for sortedBam in sortedBamList:
        filteredBam = sortedBam.replace(".bam", "_filtered.bam")

        filter_cmd = 'bamm filter -b {} --percentage_id {} --percentage_aln {} -o {}'.format(
            sortedBam, results.percent_id, results.percent_aln, destDir)
        filter_stdout, filter_stderr = execute(filter_cmd, logfile)
        logfile.write(filter_stderr)

        # BAM file from 'source' directory directory will be filtered and copied to dest dir, THAT file loc is needed
        filteredBams.append(os.path.join(destDir, os.path.basename(filteredBam)))

    return filteredBams


def samtools_sort_and_index(fileList, destDir, logfile):

    sorted_bam_files = []

    for n, files in enumerate(fileList):

        sorted_pre_bam = os.path.basename(files.replace(".bam", "_sorted.bam"))
        sorted_pre_bam_path = os.path.join(destDir, sorted_pre_bam)

        samtools_sort_cmd = 'samtools sort {} -o {}'.format(files, sorted_pre_bam_path)
        sort_stdout, sort_stderr = execute(samtools_sort_cmd, logfile)
        logfile.write('Samtools Sort')
        logfile.write(sort_stderr)

        sorted_bam = sorted_pre_bam_path
        indexed_bam = sorted_pre_bam_path + ".bai"

        samtools_index_cmd = 'samtools index {} {}'.format(sorted_pre_bam_path, indexed_bam)
        index_stdout, index_stderr = execute(samtools_index_cmd, logfile)
        logfile.write('Samtools index')
        logfile.write(index_stderr)

        sorted_bam_files.append(sorted_bam)  # Will be used later

    return sorted_bam_files


def bamm_parse(sortedBamList, logfile):

    sortbamList = " ".join(sortedBamList)

    bamm_parse_cmd = "bamm parse --coverages {} --coverage_mode {} -t {} -b {}".format(
        results.coverages, results.coverage_mode, results.threads, sortbamList)

    parse_stdout, parse_stderr = execute(bamm_parse_cmd, logfile)
    logfile.write(parse_stderr)

    return None


def bam_filter_cov(bamFileList, destDir, logfile):

    bam_cov_filtered_files = []
    for bamFile in bamFileList:

        filtered_bam_out = os.path.basename(bamFile) + '_cov-filtered.bam'
        filtered_bam_out_path = os.path.join(destDir, filtered_bam_out)
        filtered_bam_cmd = 'filter_bam_file_coverage.py --pysam_version 0.7.5 --bam {} --out {} --th {}'.format(
            bamFile, filtered_bam_out_path, results.cov_filter)

        filter_stdout, filter_stderr = execute(filtered_bam_cmd, logfile)
        logfile.write(filter_stderr)

        if 'ValueError' in filter_stderr:
            continue  # Don't write to the bam coverage file

        bam_cov_filtered_files.append(filtered_bam_out_path)

    return bam_cov_filtered_files


def normalize_coverages(coverages_fn, metagenomes_fn, logfile):

    metagenome_series = pd.read_csv(metagenomes_fn, sep=',', engine='python', index_col=0, squeeze=True)

    coverages_df = pd.read_csv(coverages_fn, sep='\t', engine='python')
    coverages_df.set_index('#contig', inplace=True)

    # Need to drop Lengths or else filtering on aligned bp will be wrong
    coverages_df.drop('Length', axis=1, inplace=True)

    # This is why length must be dropped above... also makes any potentially large calculators later less intensive
    coverages_df = coverages_df.loc[(coverages_df != 0).any(axis=1)]

    # BamM and other scripts append to names, and after a few rounds the naming can be lengthy
    # ... not totally sure this is necessary, as the nameDict (below) will find the appropriate match and rename them
    coverages_df.rename(
        columns=lambda x: os.path.basename(x).replace('_filtered_sorted_sorted.bam_cov-filtered_sorted.bam', ''),
        inplace=True)

    # The names supplied by the user through metagenomes_fn are unlikely to be identical to the read file names.... in fact,
    # they'll likely NEVER be the same (e.g., R1 and R2 versions). In response, will need to create mapping dictionary
    # Dictionary will be used to rename, and THEN can divide by the index names
    nameDict = {}

    for cov_column in coverages_df.columns:
        for sites in metagenome_series.index:

            if sites in cov_column:
                nameDict[cov_column] = sites

    coverages_df.rename(columns=nameDict, inplace=True)

    coverages_df *= 1000000000

    final_coverages_df = coverages_df.div(metagenome_series, axis='columns')

    return final_coverages_df


def create_map(dataframe):

    from palettable.colorbrewer.sequential import PuBuGn_9

    # http://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib
    cmap = PuBuGn_9.mpl_colormap
    cmap.set_under('w')
    #cmap.set_over('k')

    # Because only coverages >=75% at a site/ref is kept, it's going to have quite a bit more than 1 bp covered...
    # If using robust=True, it will ignore the extreme values... which blackens nearly all the non 0 values
    # Can't use cbar_ax, https://github.com/mwaskom/seaborn/issues/471
    ax = sns.clustermap(dataframe,
                        cmap=cmap,
                        linewidths=0.1,
                        vmin=0.01,
                        standard_scale=None,
                        col_cluster=False)  # Specialized version of heatmap, which is why ax_heatmap is called
    # vmin=0.01, vmax=500,
    # X, Y or 'both'
    ax.ax_heatmap.tick_params(axis='y', which='major', labelsize=7)
    ax.ax_heatmap.tick_params(axis='x', which='major', labelsize=7)
    ax.ax_heatmap.set_xlabel("Samples")
    ax.ax_heatmap.set_ylabel("Populations")
    #ax.ax_heatmap.set(xscale="log", yscale="log") Compresses axes to a single line

    # Set information about the colorbar, access to which is mistakenly hidden
    cax = plt.gcf().axes[-1]
    cax.set_title(label='Population relative abundance in a sample \n(Bp mapped / Kb contig / Mb metagenome)')
    cax.tick_params(size=6)

    cticks_label = [label.get_text() for label in cax.get_yticklabels()]
    cticks_values = cax.get_yticks()


    for label in ax.ax_heatmap.get_yticklabels():
        label.set_rotation(0)

    #yaxis = cax.yaxis
    #yaxis.set_ticks([0, 0.3782, 0.7602, 1.02720654])
    #cax.set_yticklabels([1, 10, 100, 500])

    for label in ax.ax_heatmap.get_xticklabels():
        label.set_rotation(45)

    return ax

if __name__ == '__main__':

    # Set up directories
    cwd = os.getcwd()

    log = open(results.log_fn, 'w')

    # Get all *.bam files
    workable_files = file_finder(results.dir, '.bam')
    log.write('Bam files found:\n')
    pprint(workable_files, log)

    # If metagenome_size is selected, need to move it within working dir context
    metagenome_size_fn = False
    if results.metagenome_size:
        metagenome_size_finder = file_finder(results.dir, '{}'.format(results.metagenome_size))
        log.write('Metagenome size file selected, copying to working directory...\n')
        metagenome_size_fn = cp_files(metagenome_size_finder, cwd, log)[0]  # If it's more than 1...

    steps = 1  # Keep track of directory creation

    sourcesDir = os.path.join(cwd, '{}.sources'.format(steps))
    if not os.path.exists(sourcesDir):
        os.makedirs(sourcesDir)
        log.write('Created directory: {}'.format(sourcesDir))
        steps += 1

    log.write('Directory content:\n')
    pprint(os.listdir(sourcesDir), log)

    # Copy files to scratch and get full file paths
    listOfBams = cp_files(workable_files, sourcesDir, log)

    log.write('Bam file locations:\n')
    pprint(listOfBams, log)

    if results.percent_id:
        log.write('Percent ID (--percent-id) specified. Will now BamM filter BAM files...\n')

        filteredBamDir = os.path.join(cwd, '{}.filtered'.format(steps))
        if not os.path.exists(filteredBamDir):
            os.makedirs(filteredBamDir)
            log.write('Created directory: {}\n'.format(filteredBamDir))
            steps += 1

        # Update BAM files
        listOfBams = bamm_filter(listOfBams, filteredBamDir, log)

        # Update BAM file location notification
        log.write('Filtered Bam Files:\n')
        pprint(listOfBams, log)

    sortedIndexedBams = os.path.join(cwd, '{}.sorted-indexed'.format(steps))
    if not os.path.exists(sortedIndexedBams):
        os.makedirs(sortedIndexedBams)
        log.write('Created directory: {}\n'.format(sortedIndexedBams))
        steps += 1

    sorted_bam_files = samtools_sort_and_index(listOfBams, sortedIndexedBams, log)

    log.write('Sorted Bam Files:\n')
    pprint(sorted_bam_files, log)

    # Simon's script to filter by sequence coverage (%), takes 1 BAM file as input, exports another user-named BAM file
    finalBams = os.path.join(cwd, '{}.coverageFiltered'.format(steps))
    if not os.path.exists(finalBams):
        os.makedirs(finalBams)
        log.write('Created directory: {}\n'.format(finalBams))
        steps += 1

    bam_cov_filter_files = bam_filter_cov(sorted_bam_files, finalBams, log)

    log.write('Coverage-Filtered Bam Files:\n')
    pprint(bam_cov_filter_files, log)

    sortedFilteredBams = os.path.join(cwd, '{}.reSorted-reFiltered'.format(steps))
    if not os.path.exists(sortedFilteredBams):
        os.makedirs(sortedFilteredBams)
        log.write('Created directory: {}\n'.format(sortedFilteredBams))
        steps += 1

    # Need to re-sort and re-index files
    bam_cov_filt_sort_files = samtools_sort_and_index(bam_cov_filter_files, sortedFilteredBams, log)

    # BamM parse command
    bamm_parse(bam_cov_filt_sort_files, log)

    # Adjust for metagenome size and coverage
    if results.metagenome_size:
        coverage_df = normalize_coverages(os.path.join(cwd, results.coverages), metagenome_size_fn, log)

        coverage_df.to_csv('adj_coverage_table.csv')

        axes = create_map(coverage_df)

        axes.savefig('adj_coverage_table.svg', dpi=600)

    pprint('Program Complete', log)