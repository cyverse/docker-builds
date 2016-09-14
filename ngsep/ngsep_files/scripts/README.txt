These scripts facilitate running the NGSEP pipeline in a command line environment. Before running these scripts, the following set up steps must be performed:

1. Create a directory called "reads" and locate your fastq files in this directory. The scripts assume that reads are paired-end, fastq files are compressed (with gzip), and for each sample the file names have the following format:
<SAMPLE_NAME>_1.fastq.gz
<SAMPLE_NAME>_2.fastq.gz

2. Create a subdirectory at the same level of the "reads" directory (I usually call it "mapping"), and place the scripts in this directory.

3. Download and decompress bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) in any directory. The scripts do not require that you install it, which is convinient for cases when you do not have root access. Edit the scripts "runTestInsertLength" and "runMapping" and change the value of the variable "BOWTIE2" to indicate the
full path where bowtie2 is located.

4. Save the reference in a single fasta file and index it with bowtie2-build using the following command:

/path/to/bowtie2/bowtie2-2.2.4/bowtie2-build reference.fa reference.fa

This can be done at any directory in the file system but this directory needs to be set changing the variable "REFERENCE" in each script. Also, create a separate text file with the chromosome (scaffold) names using the following command:

awk '{if(substr($1,1,1)==">") print substr($1,2) }' reference.fa > reference_seqNames.txt

This file will be used during the merging steps.

5. Download the jar files NGSToolsApp_<VERSION>.jar (http://sourceforge.net/projects/ngsep/files/Library/) and SortSam.jar (http://sourceforge.net/projects/picard/files/picard-tools/) and place them in one directory. Again this can be any directory but it needs to be indicated to the scripts changing the variable "JARS_DIR". The scripts "runMapping",
"runNGSEP" and "runGenotyping" do not have a version number for the file NGSToolsApp.jar. Either create a symbolic link without the version number:

ln -s NGSToolsApp_<VERSION>.jar NGSToolsApp.jar

or set the right version in the variable "NGSEP" of each script 


6. Grant execution parameters to every script using the following command;

chmod 755 run* 

Finally, the scripts assume that you are working in a 64 bit environment and that both java and bowtie2 are downloaded and run in 64-bit mode.

These are the steps we usually follow to analyze samples:

1. Estimate the insert size distribution. The script for this is "runTestInsertLength". It requires only one parameter, which is the sample name. Usage is as follows:

./runTestInsertLength <SAMPLE_NAME>

This script takes the first 250000 fragments, aligns them to the genome and provides a distribution of predicted insert lengths based on the alignments. This distribution can be
used to decide the minimum and maximum insert lengths for the next step. The script has a second optional parameter which can be used to indicate bowtie2 that base quality scores are in phred+64 format. In that case, the usage is:

./runTestInsertLength <SAMPLE_NAME> --phred64

2. Align reads to the reference and sort alignments. This scripts receives the sample name, the minimum insert length and the maximum insert length. For example, if the minimum and maximum insert lengths obtained in the previous step are 200 and 400 respectively, the command will look as follows:

./runMapping <SAMPLE_NAME> 200 400

This script also has an optional parameter to indicate bowtie2 that base quality scores are in phred+64 format. In that case, the usage is:

./runMapping <SAMPLE_NAME> 200 400 --phred64

Besides the bam files with the aligned reads and with the sorted aligned reads, this script creates the file <SAMPLE_NAME>_bowtie2_readpos.stats and <SAMPLE_NAME>_bowtie2_coverage.stats with quality statistics (differences against the reference genome) and with the coverage distribution. The quality statistics can be
inspected manually to determine how many bases should be ignored in the 5' end and/or in the 3' end respectively to call variants. This information can also be calculated using the script "calculateReadposStatsPeaks". This script calculates the average percentage of differences in the first 75 bp and then reports read positions (from 5' to 3') in which the percentage of differences with the reference is more than 2% and more than two times the average. The usage is as follows:

./calculateReadposStatsPeaks <SAMPLE_NAME>_bowtie2_readpos.stats

3. Call variants against the reference using NGSEP. The script runNGSEP has three parameters: sample name, number of bases to ignore in the 5' end and number of bases to ignore in the 3' end. For example, if the previous step suggests that we should ignore 2bp in the 5' end and 4bp in the 3' end, the command would be:

./runNGSEP <SAMPLE_NAME> 2 4

If all base calls should be taken into account, just include two zeros:

./runNGSEP <SAMPLE_NAME> 0 0

4. Merge variants from different samples. This step is used to determine the sites in which at least one sample has a genotype different than the reference. More stringent filters can be done at later stages. The command for this is MergeVariants, which is explained in the README of NGSEP and should look like this:

java -jar /path/to/jars/NGSToolsApp_<VERSION>.jar MergeVariants /path/to/reference/reference_seqNames.txt AllSamples_variants.vcf *_NGSEP.vcf

The file "AllSamples_variants.vcf" is the output file containing all genomic locations with some evidence of variation, as explained above. This file still does not contain any genotype call for any sample. The file can have any name but this name should be updated in the script "runGenotyping" (variable KNOWN_VARS) in order to run the next step.

5. Genotype samples on the sites identified in the previous step. This is done with the script runGenotyping, which has the same usage of "runNGSEP":

./runGenotyping <SAMPLE_NAME> 2 4

This should be executed sample by sample. At the end, the files *_NGSEP_gt.vcf will contain genotype calls for all sites identified in the step 4 (including reference calls).
I usually do not include quality filters at this stage to retain as much information as possible. The quality filter can be done with the FilterVCF command of NGSEP.

6. Final merge of the genotype files. This is done with the MergeVCF command of NGSEP as follows:

java -jar /path/to/jars/NGSToolsApp_<VERSION>.jar MergeVCF /path/to/reference/reference_seqNames.txt *_NGSEP_gt.vcf > AllSamples_genotypes.vcf

The file "AllSamples_genotypes.vcf" will be the final file with genotype calls for every site identified in the step 4 and for every sample analyzed. This file can be annotated, filtered and converted to other formats using NGSEP as explained in the NGSEP README (http://sourceforge.net/projects/ngsep/files/Library/)















