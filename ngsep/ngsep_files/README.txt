NGSTools - Java tools for analysis of Next Generation Sequencing (NGS) data
Version 2.1.3 (24-07-15)
===========================================================================

The NGSTools package provides an object model to enable different kinds of
analysis of Next Generation Sequencing (NGS) data. The most important
component of NGSTools is the variants detector, which performs accurate
detection and genotyping of Single Nucleotide Variants (SNVs), small and
large indels, and Copy Number Variants (CNVs). NGSTools provides utilities
to calculate quality and coverage statistics, to perform functional
annotation of variants, to merge varian calls for different samples, and
to filter and convert VCF files. The format of choice to process alignments
in every component in this package is SAM (or BAM), which allows to
integrate NGSTools with commonly used mapping programs as Bowtie2 
(http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and other analysis
packages such as SAMTools (http://samtools.sourceforge.net/) or picard
(http://picard.sourceforge.net/).

--------------------
Building NGSToolsApp 
--------------------

NGSTools has been compiled and run successfully on the standard jdk version
1.6.0. To build the distribution library NGSToolsApp.jar on a unix based
command line environment run the following commands in the directory where
NGSTools_2.1.3.tar.gz is located:

tar -xzvf NGSTools_2.1.3.tar.gz
cd NGSTools_2.1.3
make all

---------------
Asking for help
---------------

It is possible to obtain usage information for each module by typing:

java -jar NGSToolsApp_2.1.3.jar <MODULE> --help

General information and the list of modules can be obtained by typing:

java -jar NGSToolsApp_2.1.3.jar [ --help | --version | --citing ]

-------------------------------------------
Calling variants with the Variants detector
-------------------------------------------

The main module of NGSTools is the variants detector. Basic usage requires an
alignments file in SAM or BAM format, the reference genome that was used to
produce the alignments, and a prefix for the output files. The full usage is as
follows: 

USAGE: 

java -jar NGSToolsApp_2.1.3.jar FindVariants <OPTIONS> <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

OPTIONS:	
	-h FLOAT		: Heterozygosity rate. Default: 0.001
	-querySeq STRING	: Call variants just for this sequence name 
	-start INT		: Call variants just from this position in the
				  given query sequence
	-end INT		: Call variants just until this position in the
				  given query sequence
	-ignoreLowerCaseRef	: Ignore sites where the reference allele is
				  lower case. 
	-maxAlnsPerStartPos INT	: Maximum number of alignments allowed to start
				  at the same reference site. This parameter
				  helps to control false positives produced by
				  PCR amplification artifacts. Default 5 
	-s			: Consider secondary alignments while calling
				  SNVs
	-minAltCoverage	INT	: Minimum coverage of the alternative allele to
				  call a SNV. Default: 0
	-maxAltCoverage	INT	: Maximum coverage of the alternative allele to
				  call a SNV. Default: 0 (No filter)
	-minQuality	INT	: Minimum genotype quality to accept a SNV call
				  Genotype quality is calculated as 1 minus the
				  posterior probability of the genotype given
				  the reads (in phred scale). Default: 0
	-maxBaseQS INT		: Maximum value allowed for a base quality
				  score. Larger values will be equalized to
				  this value. This parameter allows to reduce
				  the effect of sequencing errors with high
				  base quality scores. Default: 0 (No limit)
	-ignore5 INT		: Ignore this many base pairs from the 5' end
				  of the reads. Default: 0
	-ignore3 INT		: Ignore this many base pairs from the 3' end
				  of the reads. Default: 0
	-ploidy INT		: Ploidy of the sample to be analyzed. 
				  Default 2
	-knownSVs FILE		: File with coordinates of known structural
				  variants in GFF format.
	-genomeSize INT		: Total size of the genome to use during
				  detection of CNVs. This should be used when
				  the reference file only includes a part of
				  the genome (e.g. a chromosome or a partial
				  assembly)   
	-binSize INT		: Size of the bins to analyze read depth.
				  Default: 100
	-algCNV	STRING		: Comma-separated list of read depth algorithms
				  to run. Default: CNVnator
	-maxPCTOverlapCNVs INT	: Maximum percentage of overlap of a new CNV
				  with an input CNV to include it in the output
				  Default: 100 (No filter)
	-maxLenDeletion INT	: Maximum length of deletions that the read-pair
				  analysis can identify. Default: 1000000
	-ignoreProperPairFlag	: With this option, the proper pair flag will
				  not be taken into accout to decide if the ends
				  of each fragment are properly aligned. By
				  default, the distribution of insert length is
				  estimated only taking into account reads with
				  the proper pair flag turned on
	-minSVQuality INT	: Minimum quality score (in PHRED scale) for
				  structural variants. Default: 20
	-sizeSRSeed INT		: Size of the seed to look for split-read
				  alignments. Default: 8
	-ignoreXS		: Ignores the optional field XS to decide if an
				  alignment is unique. While bowtie2 only
				  outputs the XS field for reads with multiple
				  alignments, BWA-MEM and BWA-SW outputs this
				  field for every alignment in the BAM file.
				  Hence, this flag should always be used to
				  process BAM files made with BWA-MEM and 
				  BWA-SW
	-sampleId STRING	: Id of the sample that will appear in the
				  output vcf file 
	-psp			: Flag to print a header in the VCF file with
				  the id and the ploidy of the sample			  
	-knownVariants FILE	: VCF file with variants to be genotyped. Only
				  these variants will appear in the output vcf
				  file. With this option homozygous calls to
				  the reference allele will be reported  
	-genotypeAll		: Report all covered sites in the genome.
				  CAUTION: For WGS samples, and even for 
				  RAD/GBS samples in organisms with medium to
				  large genomes this produces a huge file. 
	-noRep			: Turn off finding repetitive regions based on
				  reads with multiple alignments
	-noRD			: Turn off read depth analysis to identify CNVs
	-noRP			: Turn off read pair analysis to identify large
				  indels and inversions
	-noSNVS			: Turn off SNV detection
	-noNewCNV		: Turn off finding new CNVs with the read depth
				  analysis. Input CNVs and repeats will still
				  be genotyped using the RD distribution

Alignments should be provided in SAM V-0.1.2 format
(see http://samtools.sourceforge.net for details) or its binary version (BAM).
The reference genome should be provided in fasta format. It is assumed that
the sequence names in the alignments file correspond with the sequence names in
this reference assembly. The output for SNVs and small indels is a VCF file
(see http://www.1000genomes.org/node/101 for details). VCF files produced by
the variants detector include the custom field CNV in the INFO column for each
SNV or small indel, which in general reports the number of samples in which a
CNV covers the SNV or small indel. SNV calls also include the genotype field
AAC which will report the number of As, Cs, Gs , and Ts observed in the reads.

Structural variants are reported in a gff file
(see http://www.sequenceontology.org/gff3.shtml for details). Besides the
genomic coordinates and the quality information, this gff will contain the
algorithm that originated each variant call, which allows to differentiate
calls from the repetas detector, the read depth analysis, and the read pair
analysis. This gff can be directly loaded as a custom annotation file in
gbrowse or in the UCSC genome browser to visualize the genomic regions with
structural variants and to relate them with functional elements. This file
can also be used as a parameter of the variants detector (option "-knownSVs")
which is useful to avoid recalculation of structural variants while genotyping
known variants. 

WARNING: Default parameters of the variants detector are designed to maximize
the number of called variants which will generate false positives in samples
with small coverage or high error rates. For conservative SNV calling use:

java -jar NGSToolsApp_2.1.3.jar FindVariants -ignoreLowerCaseRef -maxAlnsPerStartPos 2 -minQuality 40 -maxBaseQS 30 <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

If the error rate towards the three prime end increases over 2% you can also
use the option -ignore3 to ignore errors at those read positions or use the
Clip command (see below)

WARNING 2: For RAD Sequencing or GBS samples, using the default value of the
parameter to control for PCR duplicates (maxAlnsPerStartPos) will yield very
low sensitivity. We recommend to increase the value of the parameter to about
100 to retain high sensitivity while avoiding a severe penalty in memory usage.
Also, structural variants should not be called using these data. The usage for
conservative variant calling in RAD-Seq or GBS samples becomes:

java -jar NGSToolsApp_2.1.3.jar FindVariants -noRep -noRD -noRP -ignoreLowerCaseRef -maxAlnsPerStartPos 100 -minQuality 40 -maxBaseQS 30 <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>
  
---------------------------------
Calculating mismatches statistics
---------------------------------

This module takes one or more sets of alignments and a reference genome and
counts the number of mismatches with the reference for each read position from
5' to 3' end. This report is useful to detect sequencing error biases. The
usage for this tool is the following

USAGE:

java -jar NGSToolsApp_2.1.3.jar QualStats <OPTIONS> <REFERENCE_FILE> <ALIGNMENTS_FILE>*

OPTIONS:
	-ignoreXS	: Ignores the optional field XS to decide if an
			  alignment is unique. See the same option in
			  the variants detector for more details
						  
The file(s) with alignments must be given in SAM or BAM format and the
reference file in fasta format. The output is a text file with five columns:
- Position: 1- based from 5' to 3'
- Number of reads with a base call different than the reference (Considering
  all alignments)
- Number of reads with a base call different than the reference (Considering
  only reads with unique alignments)
- Number of total alignments counted with read length equal or larger than the
  position in the first column. The percentage of mismatches including all
  alignments is the ratio of column 2 divided by this column
- Number of uniquely aligned reads counted with read length equal or larger
  than the position in the first column. The percentage of mismatches for
  uniquely aligned reads is the ratio of column 3 divided by this column

-------------------------------
Calculating coverage statistics
-------------------------------

This module calculates the number of base pairs that are covered by reads at
each coverage level from 1 to a maximum. This statistic is useful to visualize
how uniform was the sequencing process over the genome. The usage is as follows

USAGE:

java -jar NGSToolsApp_2.1.3.jar CoverageStats <ALIGNMENTS_FILE> <OUTPUT_FILE>


The alignments file must be given in SAM or BAM format. The output is a text
file with three columns:
- Coverage
- Number of reference sites with this coverage (Considering all alignments)
- Number of reference sites with this coverage (Considering only reads with 
  unique alignments)

---------------------------------
Functional annotation of variants
---------------------------------

This module takes a VCF file produced by NGSEP, the reference genome in fasta
format, and a gff3 file with gene annotations related with the given genome
(see http://www.sequenceontology.org/gff3.shtml for details) and generates a
VCF file which includes the functional information related with each variant.
The usage is as follows

USAGE:

java -jar NGSToolsApp_2.1.3.jar Annotate <OPTIONS> <VARIANTS_FILE> <TRANSCRIPTOME_MAP> <REFERENCE_FILE>

OPTIONS:
	-u INT	: Maximum bp before a gene to classify a variant as Upstream.
		  Default: 1000
	-d INT	: Maximum bp after a gene to classify a variant as Downstream.
		  Default: 300

The vcf file with functional annotations is written in the standard output.
Annotations are included using the following custom fields in the INFO column:

TA (STRING): 	Annotation based on a gene model. Possible annotations include
		Intergenic, Intron, FivePrimeUTR, ThreePrimeUTR, Upstream, 
		Downstream, NCRNA, Synonymous, Missense, Nonsense, Frameshift,
		and ExonJunction
TID (STRING):	Id of the transcript related with the gene annotation in the TA
		field
TGN (STRING): 	Name of the gene related with the annotation in the TA field
TCO (FLOAT):	For variants in coding regions, position in the aminoacid
		sequence where the variant is located. The integer part is the
		1-based position of the mutated codon. The decimal part is the
		codon position.

------------------------
Clipping read alignments
------------------------

NGSTools includes a utility to clip a given number of base pairs from the 5'
end and from the 3' end of each alignment in a SAM file. The input for this
tool is a file of read alignments in SAM or BAM format, the number of bases to
clip from the 5' end and the number of bases to clip from the 3' end. The usage
is the following:

USAGE: 

java -jar NGSToolsApp_2.1.3.jar Clip <OPTIONS> <ALIGNMENTS_FILE> <OUTPUT_FILE> <CLIPPING_5PRIME_END> <CLIPPING_3PRIME_END>


The only option is -s which lets the program know that alignments are sorted by
reference position. In this mode, output alignments will also be sorted by
reference position 

----------------------------------------
Merging variants from individual samples
----------------------------------------

NGSTools can be used to merge variants from different samples into an
integrated VCF file. The pipeline for this purpose is as follows.

The first step is to generate a file including the whole set of variants called
in at least one of the samples. This can be done calling the MergeVariants
command as follows:

USAGE:

java -jar NGSToolsApp_2.1.3.jar MergeVariants <SEQUENCE_NAMES_FILE> <OUTPUT_FILE> <VARIANTS_FILE>*


The sequence names file is a text file which just has the ids of the sequences
in the reference. It is used by the program to determine the order of the
reference sequences. In unix systems this file can be obtained running the
following command on the fasta file with the reference genome:

awk '{if(substr($1,1,1)==">") print substr($1,2) }' <REFERENCE_FILE> > <SEQUENCE_NAMES_FILE>

The output file of this merge program is a vcf with the union of variants
reported by the input files but without any genotype information. 

The second step is to genotype for each sample the variants produced at the
first step using the variants detector (See FindVariants command). For each
sample, the command to execute at this stage (in conservative mode and assuming
WGS data) should look like this:

java -jar NGSToolsApp_2.1.3.jar FindVariants -noRep -noRD -noRP -ignoreLowerCaseRef -maxAlnsPerStartPos 2 -minQuality 40 -maxBaseQS 30 -knownSVs <SVS_FILE> -knownVariants <VARS_FILE> <REFERENCE> <INPUT_FILE> <OUTPUT_PREFIX>

where SVS_FILE is the file with structural variation for the sample obtained
during the first run of the variants detector, and VARS_FILE is the output file
obtained in the first step of the merging process. At the end, this will
produce a second set of vcf files which will differ from the first set in the
sense that they will include calls to the reference allele. The third step is
to join these new vcf files using the following command:

USAGE:

java -jar NGSToolsApp_2.1.3.jar MergeVCF <SEQUENCE_NAMES_FILE> <GENOTYPED_VARIANTS_FILE>*


This command will write to standard output the final vcf file with the genotype
calls for each variant on each sample.

For organisms with small genomes, a VCF for the whole genome can be generated
for each sample using the option -genotypeAll in the variants detector. Then,
only the third step will be needed to mix the VCF files generated in this mode.

-------------------
Filtering VCF files
-------------------

This module implements different filters on VCF files with genotype
information. It writes to standard output a VCF file with variants passing the
filtering criteria. Since version 2.0.6, the default behavior does not perform
any filtering. The filtering order is as follows: first, it executes the
distance filter (-d option), then the filtering of samples and genotypes (-saf,
-fs, -q and -minC options). Finally, it recalculates the number of samples 
genotyped, the number of alleles called and the MAF to execute the remaining 
filters.

USAGE:

java -jar NGSToolsApp_2.1.3.jar FilterVCF <OPTIONS> <INPUT_FILE>

OPTIONS:
	-frs FILE	: File with genomic regions in which variants should be
			  filtered out. The format of this file should contain
			  at least three columns: Sequence name (chromosome),
			  first position in the sequence, and last position in
			  the sequence. Both positions are assumed to be
			  1-based.
	-srs FILE	: File with genomic regions in which variants should be
			  selected. The format of this file should contain at
			  least three columns: Sequence name (chromosome),
			  first position in the sequence, and last position in
			  the sequence. Both positions are assumed to be
			  1-based.
	-d INT		: Minimum distance between variants.
	-g FILE		: File with the reference genome to calculate the
			  GC-Content of the region surrounding the variant.
	-minGC FLOAT	: Minimum percentage of GC of the 100bp region
			  surrounding the variant.
	-maxGC FLOAT	: Maximum percentage of GC of the 100bp region
			  surrounding the variant.
	-q INT		: Minimum genotyping quality score (GQ field for each
			  genotype call in the vcf file).
	-s		: Keep only biallelic SNVs.
	-fi		: Flag to filter sites in which only one allele is
			  observed.
	-fir		: Flag to filter sites in which only the reference
			  allele is observed.
	-fia		: Flag to filter sites in which only one alternative
			  allele is observed.
	-minI INT	: Minimum number of individuals genotyped to keep the
			  variant.
	-minC INT	: Minimum coverage to keep a genotype call.
	-minMAF FLOAT	: Minimum minor allele frequency over the samples in
			  the VCF file.
	-maxMAF FLOAT	: Maximum minor allele frequency over the samples in
			  the VCF file.
	-maxCNVs INT	: Maximum number of samples with copy number variation
			  in the region where the variant is located.
	-gene STRING	: Id of the gene or the transcript related with the
			  variant.
	-a STRING	: Types of functional annotations (Missense, Nonsense,
			  Synonymous, etc) related with the variant. More than
			  one annotation can be set as a comma-separated list
	-saf FILE	: File with the ids of the samples to be selected (or
			  filtered, see -fs option). The file should have one
			  line per sample, being the first column the sample
			  id. Other columns in the file are ignored.
	-fs		: Flag to filter the samples provided with the -saf
			  option instead of selecting them. 

----------------------------------
Convert VCF files to other formats
----------------------------------

This module allows to convert genotype calls in VCF format to other formats
commonly used to perform different kinds of analysis.

USAGE:

java -jar NGSToolsApp_2.1.3.jar ConvertVCF <OPTIONS> <INPUT_FILE> <OUTPUT_PREFIX>

OPTIONS:
	-printStructure		: Prints input format for structure
	-printFasta		: Prints a virtual multiple sequence alignment
				  in fasta format. Useful to build phylogenetic
				  trees
	-printrrBLUP		: Prints the input files for rrBLUP
	-printMatrix		: Prints genotypes in a simple ACGT format
				  which can be imported to excel 
	-printHapmap		: Prints Hapmap format, which can be used in
				  programs such as Tassel
	-printSpagedi		: Prints the input files for Spagedi
	-printPlink		: Prints the input files for Plink
	-printHaploview		: Prints the input files for Haploview
	-printEmma		: Prints the input files for Emma
	-printPowerMarker	: Prints the input files for Powermarker
	-printEigensoft		: Prints the input files for Eigensoft
	-printFlapjack		: Prints the input files for Flapjack
	-printDarwin		: Prints the input files for DarWin
	-printTreeMix		: Prints the input files for TreeMix
	-printJoinMap		: Prints the input file to build genetic maps
				  with JoinMap
	-printPhase		: Prints the input file for PHASE
	-p1 STRING		: Id of the first parent for conversion to
				  JoinMap
	-p2 STRING		: Id of the second parent for conversion to
				  JoinMap
	-s STRING		: Name of the sequence (chromosome) for
				  conversion to PHASE 
	-p FILE			: File with population assignments for the
				  samples. This should be a two column text
				  file with the sample ids in the first column
				  and the ids of the populations in the second
				  column. Required for conversion to TreeMix

WARNING: FASTA convertion does not use IUPAC codes, heterozygous SNPs are 
changed to N.
WARNING 2: Plink is only designed for humans, therefore it will only work for
22 sequences (chromosomes). If a sample exceeds this number, it is convenient 
to reduce the number of chromosomes and to remove all scaffolds.
WARNING 3: To generate dendograms in Tassel, it is better to use the HapMap 
format.

------------------------------
Calculating summary statistics
------------------------------

This module writes to the standard output a report with the numbers of variants
included in a VCF file for different categories. Although it can be called for
any VCF file generated by the pipeline, this report is specially useful when a
complete population is being processed and merged into a single annotated file.

USAGE:

java -jar NGSToolsApp_2.1.3.jar SummaryStats <OPTIONS> <INPUT_FILE> 

OPTIONS:
	-m INT			: Minimum number of samples genotyped to
				  accurately calculate the minor allele
				  frequency. Default: 20

-----------------------------------------
Calculating diversity statistics per site
-----------------------------------------

This module produces basic diversity statistics for each variant in a VCF file.
It receives a VCF file and an optional text file with population assignments for
each sample included in the VCF and prints to the standard output the
coordinates of each variant plus the following statistics separated by
semicolon:

1. Number of samples genotyped
2. Expected heterozygosity (under HWE)
3. Observed heterozygosity
4. F-statistic (1-OH/EH)
5. Minor allele frequency (MAF)
6. Chi-square value of departure from HWE
7. Uncorrected p-value of the Chi-square test for departure from HWE

If a file with population assignments is provided, this module will output one
column of statistics for the whole group and one column for each population.

USAGE:

java -jar NGSToolsApp_2.1.3.jar DiversityStats <INPUT_FILE> <POPULATIONS_FILE>


The populations file is a tab-delimited text file with two columns: sample id
and population id.

-------------------
Comparing VCF files
-------------------
This module allows to compare the genotype calls included in two different VCF
files, calculating the number and percentage of homozygous and heterozygous
diffrences between every pair of samples. It write to standard output a 
tab-delimited report with the following fields:

1. Id sample VCF 1
2. Id sample VCF 2
3. Number of variants genotyped in sample 1
4. Number of variants genotyped in sample 2
5. Number of variants genotyped in both samples
6. Number of heterozygous differences
7. Percentage of heterozygous differences (sixth field / fifth field)
8. Number of homozygous differences
9. Percentage of homozygous differences (eighth field / fifth field)
10. Number of total differences
11. Percentage of total differences (tenth field / fifth field) 

USAGE:

java -jar NGSToolsApp_2.1.3.jar CompareVCF <OPTIONS> <SEQUENCE_NAMES_FILE> <FIRST_VCF_FILE> <SECOND_VCF_FILE>

OPTIONS: 

	-g FLOAT	: Minimum percentage (0-100) of variants genotyped in
			  both samples. Default: 50.
	-d FLOAT	: Maximum percentage (0-100) of differences between the
			  pair of samples. Default: 5.

---------------------------------
Genotype imputation (In progress)
---------------------------------
This module allows imputation of missing genotypes from unphased multilocus
SNP genotype data in a VCF. The current version is a reimplementation of the
Hidden Markov Model (HMM) implemented in the package fastPHASE
(http://stephenslab.uchicago.edu/software.html). This implementation allows to
process VCF files and produces its output also as a VCF. Because the current
implementation treats the input genotypes in the vcf as haploid data, this
module is only able to impute genotypes in populations of samples with low
heterozygosity (autogamous or double-haploids). Only biallelic SNPs are imputed
and included in the output VCF file. The current implementation of the model
does not produce heterozygous imputed genotypes. We expect to include more
general cases and improve the imputation accuracy in future versions.

USAGE:

java -jar NGSToolsApp_2.1.3.jar ImputeVCF <OPTIONS> <VCF_FILE> <OUT_PREFIX>

OPTIONS: 

	-p STRING	: Comma-separated list of sample ids of the parents of
			  the breeding population. This should only be used for
			  bi-parental or multi-parental breeding populations.
	-k INT		: Maximum number of groups in which local haplotypes
			  will be clustered. See (PMID:16532393) for details
			  of the HMM implemented in the fastPHASE algorithm.
			  For bi-parental or multi-parental breeding
			  populations please set explicitly the number of
			  parents of the population even if the list of parents
			  is provided with the -p option. This allows to take
			  into account cases of populations in which some of
			  the parents are missing. Default: 20 
	-c FLOAT	: Estimated average number of centiMorgans per Kbp on
			  euchromatic regions of the genome. This value is used
			  by the model to estimate initial transitions between
			  the states of the HMM. Typical values of this
			  parameter are 0.001 for human populations, 0.004 for
			  rice and 0.35 for yeast populations. We expect to
			  implement an option to allow setting the estimated
			  recombination rate per site in future versions.
			  Default: 0.001
	-t		: If set, transition probabilities in the HMM will NOT
			  be updated during the Baum-Welch training of the HMM. 
			  Not recommended unless the -c option is set to a
			  value allowing a reasonable initial estimation of the
			  transition probabilities.


This module outputs two files, the first is a VCF file including the imputed
genotypes for the datapoints having an undecided genotype call in the input
file. The second outputs for each SNP and each sample the index of the parent
that most likely originated the observed haplotype of the individual.

-------------------
Deconvoluting reads
-------------------

This option allows to build individual fastq files for different samples from
a single file containing the reads for a whole sequencing lane in which several
samples were barcoded and sequenced. Up to this point only single read data can
be deconvoluted.

USAGE:

java -jar NGSToolsApp_2.1.3.jar Deconvolute <OPTIONS> <INDEX_FILE> <FLOWCELL> <LANE> <INPUT_FILE>

OPTIONS: 
	-o DIRECTORY	: Directory where the output fastq files will be saved
	-t STRING	: If this sequence is found within a read, the read
			  will be trimmed up to the start of this sequence
	-u		: Output uncompressed files

INDEX_FILE is a tab-delimited text file with four columns: flowcell, lane,
barcode and sampleID. Although the same index file can be used for different
lanes, each process will only deconvolute reads from the lane determined by the
required parameters FLOWCELL and LANE. If INPUT_FILE is not present or is equal
to "-", the lane reads will be read from the standard input. The zcat command
can be used then to avoid uncompressing the input file:

zcat <INPUT_FILE> | java -jar NGSToolsApp_2.1.3.jar Deconvolute <OPTIONS> <INDEX_FILE> <FLOWCELL> <LANE>

------------------------------------
Comparing read depth between samples
------------------------------------

This function compares the read depth of two samples. It takes two alignment
files and a reference genome, splits the genome into windows, and for each
window compares the read depth between the two samples. It outputs a text file
containing the list of windows of the genome in which the normalized read depth 
ratio between the two samples is significantly different from 1. The text file
contains the following columns:

1. Chromosome
2. Window start
3. Window end
4. Read depth sample 1
5. Read depth sample 2
6. Normalized read depth ratio
7. P-value

USAGE:

java -jar NGSToolsApp_2.1.3.jar CompareRD <OPTIONS> <ALIGNMENTS_FILE_1> <ALIGNMENTS_FILE_2> <REFERENCE> <OUT_PREFIX>

OPTIONS:

	-binSize INT	: Window size to be used during the read depth
			  comparison. Default: 100
	-p FLOAT	: Maximum p-value. Only the windows with a p-value
			  lower than that specified will be reported.
			  Default: 0.001
	-w		: Output an entry for every window in the genome
	-g		: Perform GC-correction of the read depth
	-b		: Perform the Bonferroni correction for multiple 
			  testing

------------------------------
Citing and supporting packages
------------------------------

A manuscript with the description of the latest modules of NGSEP has been
recently pulbished in Nuceic Acids research:

Duitama J, Quintero JC, Cruz DF, Quintero C, Hubmann G, Foulquie-Moreno MR, Verstrepen KJ, Thevelein JM, and Tohme J. (2014). 
An integrated framework for discovery and genotyping of genomic variants from high-throughput sequencing experiments. 
Nucleic Acids Research. doi: 10.1093/nar/gkt1381
http://nar.oxfordjournals.org/content/early/2014/01/11/nar.gkt1381.full

Details of variant detection algorithms implemented in NGSEP can be found in
the following publications:

SNV detection:
Duitama J, Srivastava PK, and Mandoiu II. (2012). 
Towards accurate detection and genotyping of expressed variants from whole transcriptome sequencing data. 
BMC Genomics, 13(Suppl 2), S6. doi:10.1186/1471-2164-13-S2-S6

CNV detection (Read depth analysis):
Abyzov, A., Urban, A. E., Snyder, M., and Gerstein, M. (2011). 
CNVnator: an approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing.
Genome research, 21(6), 974–84. doi:10.1101/gr.114876.110

Yoon S, Xuan Z, Makarov V, Ye K, Sebat J. (2009).
Sensitive and accurate detection of copy number variants using read depth of coverage. 
Genome Research Sep;19(9):1586-92. Epub 2009 Aug 5. PubMed PMID: 19657104; PubMed Central PMCID: PMC2752127.

Genotype imputation:
Scheet, P and Stephens, M. (2006).
A Fast and Flexible Statistical Model for Large-Scale Population Genotype Data: Applications to Inferring Missing Genotypes and Haplotypic Phase.
American Journal of Human Genetics 78: 629-644.

Read Depth comparison:
Xie C and Tammi MT. (2009).
CNV-seq, a new method to detect copy number variation using high-throughput sequencing.
BMC Bioinformatics 10:80.

NOTE: Since version 2.1.2, we implemented a new model to integrate paired-end and split-read analysis for detection of large indels.
Benchmarking with other tools is in progress.

NGSEP is also supported by the following open source software packages:

Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
Picard: http://picard.sourceforge.net/
Jsci: http://jsci.sourceforge.net/
XChart: http://xeiam.com/xchart/
