# Docker-BowtieBatch

## Overview

**Docker-BowtieBatch** is a python-wrapper for [Bowtie2](bowtie-bio.sourceforge.net/bowtie2/), designed to run within a Docker container. The output produced by this script is suitable for use with [Docker-Read2RefMapper](https://bitbucket.org/bolduc/docker-read2refmapper).

## Installation

Since everything runs within a Docker container, there is no need to install any dependencies, other than [Docker](https://www.docker.com/).

### Using the Dockerfile

Build the Docker image from the Dockerfile:

```
docker build -t bbolduc/docker-bowtiebatch:dev .
```

The "." in the code above means to use the Dockerfile from the current directory. (the file is literally named "Dockerfile")

Now you can "run" the Docker image - almost as if it's another program on your system (typical run command).

```
docker run --rm -v /path/to/fasta/files/dir:/inputDir -w /inputDir bbolduc/docker-bowtiebatch:dev --reads-dir /inputDir
```

Note: Commands after the "bbolduc/docker-bowtiebatch:dev" are being passed to the Dockerized python script.

### Inputs

Recent improvement requests included the ability to use both paired and unpaired reads. This has now been done with the "mixed" option. For this to work, this script naturally sorts files according to their names and uses the Levenshtein distance between filenames, essentially operating under the assumption that unpaired reads follow paired reads in their naming, i.e.:

```
Reads1_F_paired.fastq
Reads1_F_unpaired.fastq
Reads1_R_paired.fastq
Reads1_R_unpaired.fastq
```

Given the above example ("mixed" option and "dist" of 3), the script would group the Reads1 into:

```
Paired:	[Reads1_F_paired.fastq
		 Reads1_R_paired.fastq]
Unpaired: [Reads1_F_unpaired.fastq
		   Reads1_R_unpaired.fastq]
```
Notice that the "distance" between the 2 paired files is 1 - there's only 1 letter difference (the F and R). For the unpaired, this is also true (the F and R). *However* the distance between the paired and unpaired *of the same direction* is 2 ("un" in unpaired). So to group all the Reads1 you need a distance of 3; +1 to account for F/R and +2 for paired/unpaired.

However, it's recognized that people (and sequencing centers) name their files differently, so users can adjust the distance to ensure correct grouping. However, a fundamental assumption for this entire process is that **filenames within a read group must be more similar within the group than between groups.** If reads don't correctly group you may have to 1) re-name files or 2) run the script twice for paired and unpaired.

File inputs can be fasta or fastq-formatted, and can be uncompressed or compressed as ".tar.gz" or ".gz" *No other options are currently allowed.*

### Minimal command line options

**--rm**: Remove the Docker container after running.

**-v**: Bind directory from the *host* system to the container system. This is the **absolute host path**, followed by **:**, then where (on the *container's* system) to mount it.

**-w**: Sets working directory *within the container*.

**--reads-dir** / **--reads**: Directory containing reads file(s) to be aligned. For ease of use, this **must be** the same directory that is mounted under the container.

**--input-db** / **--db** / **-d**: NAME of FASTA-formatted sequence file containing sequences/references to be aligned against. This must be a single file.

**--input-format** / **--fmt** / **-f**: File format of reads to be aligned. Compression is recognized (must be *.gz or *.tar.gz). Must be fasta, fastq or fq.

**--interleaved** / **-s**: Will treat *each* and *all* files as interleaved, where both forward and reverse reads are present in the same file, alternating.

### General Options

**--read-types**: Whether or not reads are paired, unpaired or mixed. See "Inputs" above for more information. Default: mixed.

**--dist**: Levenshtein distance between filenames. For paired reads, where only 1 letter of difference separates filenames, use 1. For paired/unpaired where distances could be 2-3, use 2-3. Default: 1. 

**--exclude** / **-x**: A list of strings (with -x for each item) to filter out files within the input directory. Useful for excluding other files within the directory. Not currently enabled on Stampede. Default: .py.

**--unpair-terms** / **-u**: A list of strings (with -u for each item) to help the script identify unpaired reads. By default, unpaired files are identified by their natural sort order against (assumed to be present) paired-end files.

**--pair-terms** / **-p**: A list of strings (with -p for each item) to help the script identify paired reads. By default, paired files are identified by their natural sort order against (assumed to be present) unpaired files.

**--keep-sam**: SAM files will be preserved during BAM file generation. Unless enabled, no SAM files will be generated. Default: off.

**--merge-results**: Output from all runs will be combined into a single file.

**--merge-name**: Filename to use for merged results. Not used if --merge-results is not enabled.

### Bowtie2 Options

**--alignment_type**: Whether reads should be aligned end-to-end (global) or local. Default: end-to-end

**--end-to-end_presets**: Presets for end-to-end alignments. Please see the Bowtie2 manual for parameters.

**--local_presets**: Presets for local alignments. Please see the Bowtie2 manual for parameters.

**--non-deterministic**: Flag to enable Bowtie2's using current time to re-initialize the pseudo-random number generator. Useful when the input consists of many identical reads with read names very similar to each other. Default: off.

**--trim5**: Trim X bases from 5'/left end of reads. Default: 0.

**--trim3**: Trim X bases from 3'/right end of reads. Default: 0.

**--minins**: Minimum fragment length. Default: 0.

**--maxins**: Maximum fragment length. Default: 2000.


### Authors

* Ben Bolduc

For a full list of contributions, see the Contributors file.