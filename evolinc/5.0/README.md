# EVOLINC-I 1.0
Evolinc-I is a long intergenic noncoding RNA (lincRNA) identification workflow that also facilitates genome browser visualization of identified lincRNAs and downstream differential gene expression analysis. 

Evolinc-I minimally requires the following input data

1. A set of assembled and merged transcripts from Cuffmerge or Cuffcompare in gene transfer format (GTF)
2. A reference genome (FASTA)
3. A reference genome annotation (GFF)

Optional input data

1. Transposable Elements database (FASTA)
2. Known LincRNA (GFF)
3. Transcription start site coordinates (BED)
 
# Availablility
### Using Docker image

Since there are several dependencies (these can be seen in [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-i/~/dockerfile/)) for running Evolinc-I on your linux or MAC OS, we highly recommend using the available Docker image for [Evolinc-I](https://hub.docker.com/r/cyverse/evolinc-i/) or the [Dockerfile](https://hub.docker.com/r/cyverse/evolinc-i/~/dockerfile/) to build an image and then use the built image.

```
# Pull the image from CyVerse Dockerhub
docker pull cyverse/evolinc-i:1.0

# See the command line help for the image
docker run cyverse/evolinc-i:1.0 -h 

# Run Evolinc-I on the test data. The sample data can be found in the sample_data folder in this repo
docker run --rm -v $(pwd):/working-dir -w /working-dir cyverse/evolinc-i:1.0 -c Sample_cuffcompare_out.gtf -g TAIR10_chr1.fasta -r TAIR10_chr1_genes.gff -o test_out -n 4
```

### Using CyVerse Discovery Environment

The [Evolinc-I app](https://de.cyverse.org/de/?type=apps&app-id=e980754e-8050-11e6-97c3-008cfa5ae621&system-id=de) is currently integrated in CyVerse’s Discovery Environment (DE) and is free to use by researchers. The complete tutorial is available at this [CyVerse wiki](https://wiki.cyverse.org/wiki/display/TUT/Evolinc+in+the+Discovery+Environment). CyVerse's DE is a free and easy to use GUI that simplifies many aspects of transcriptome assembly. If you do not currently have access to a high performance computing cluster, consider taking advantange of the DE.

# Issues
If you experience any issues with running Evolinc-I (DE app or source code or Docker image), please open an issue on this github repo.

# Copyright free
The sources in this [Github](https://github.com/Evolinc/Evolinc-I) repository, are copyright free. Thus you are allowed to use these sources in which ever way you like. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license.

# Citing Evolinc-I
Evolinc-I manuscript is currently under review but is available as a [bioRxiv](http://biorxiv.org/content/early/2017/02/20/110148) preprint.
