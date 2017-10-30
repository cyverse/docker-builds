# EVOLINC-I: A rapid Long-Intergenic Noncoding RNA (lincRNA) detection pipeline

# Introduction

Evolinc-I is a long intergenic noncoding RNA (lincRNA) identification workflow that also facilitates genome browser visualization of identified lincRNAs and downstream differential gene expression analysis. 

Evolinc-I minimally requires the following input data

1. A set of assembled and merged transcripts from Cuffmerge or Cuffcompare in gene transfer format (GTF)
2. A reference genome (FASTA)
3. A reference genome annotation (GFF/GTF/GFF3)

Optional input data

1. Transposable Elements database (FASTA)
2. Known LincRNA (GFF)
3. Transcription start site coordinates (BED)
 

# Availablility
### Using Docker image

Since there are several dependencies (these can be seen in [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-i/~/dockerfile/)) for running Evolinc-I on your linux or MAC OS, we highly recommend using the available Docker image for [Evolinc-I](https://hub.docker.com/r/cyverse/evolinc-i/) or the [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-i/~/dockerfile/) to build an image and then use the built image. Docker can be installed on any of three platform using the instructions from [Docker website](https://docs.docker.com/engine/installation/). You can also try [Play-With-Docker](http://labs.play-with-docker.com/) for running Evolinc-I using the below instructions 

```
# Pull the image from CyVerse Dockerhub
docker pull evolinc/evolinc-i:1.5.1
```

```
# See the command line help for the Docker image
docker run evolinc/evolinc-i:1.5.1 -h 
```

```
# Download some sample data 
git clone https://github.com/Evolinc/Evolinc-I.git
cd Evolinc-I
```

```
# Run Evolinc-I With mandatory files
docker run --rm -v $(pwd):/working-dir -w /working-dir evolinc/evolinc-i:1.5.1 -c Sample_cuffcompare_out.gtf -g TAIR10_chr1.fasta -r TAIR10_chr1_genes.gff -o test_out -n 4
```

```
# Run Evolinc-I With both mandatory and optional files
docker run --rm -v $(pwd):/working-dir -w /working-dir evolinc/evolinc-i:1.5.1 -c Sample_cuffcompare_out.gtf -g TAIR10_chr1.fasta -r TAIR10_chr1_genes.gff -b TE_RNA_transcripts.fa -t Sample_TSS_data.gff -x Sample_known_lncRNAs.gff -o test_out -n 4
```

### Using CyVerse Discovery Environment

The [Evolinc-I app](https://de.cyverse.org/de/?type=apps&app-id=e980754e-8050-11e6-97c3-008cfa5ae621&system-id=de) is currently integrated in CyVerse’s Discovery Environment (DE) and is free to use by researchers. The complete tutorial is available at this [CyVerse wiki](https://wiki.cyverse.org/wiki/display/TUT/Evolinc+in+the+Discovery+Environment). CyVerse's DE is a free and easy to use GUI that simplifies many aspects of running bioinformatics analyses. If you do not currently have access to a high performance computing cluster, consider taking advantange of the DE.

#### Step-by-step walkthroughs

Step-by-step walkthrough for running Evolinc-I on DE is available [here](https://drive.google.com/open?id=0B-ferWixi_V3cmh0QzhJeXRXSE0).
Step-by-step walkthrough for the command-line, with directions on how to change parameters is coming soon. Information on how to easily create a Cuffmerge/Cuffcompare input file from 1-many SRA IDs within the DE can be found [here](https://drive.google.com/open?id=0B-ferWixi_V3NjVpdENLUXhLZjQ)


# Issues
If you experience any issues with running Evolinc-I (DE app or source code or Docker image), please open an issue on this github repo.

# Copyright free
The sources in this [Github](https://github.com/Evolinc/Evolinc-I) repository, are copyright free. Thus you are allowed to use these sources in which ever way you like. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license.

# Citing Evolinc-I
If you have used Evolinc-I manuscript in your research, please cite as below..

*Andrew D. Nelson&ast;, Upendra K. Devisetty&ast;, Kyle Palos, Asher K. Haug-Baltzell, Eric Lyons, Mark A. Beilstein (2017). "Evolinc: a comparative transcriptomics and genomics pipeline for quickly identifying sequence conserved lincRNAs for functional analysis". Frontiers in Genetics. 1(10)*
