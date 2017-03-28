#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


my (@file_query, $database_path, $user_database_path, $annotation_path,
$user_annotation_path, @file_query2, $file_type, $STAR_ARGS);


GetOptions( "file_query=s"      => \@file_query,
            "file_query2=s"     => \@file_query2,
            "user_database=s"   => \$user_database_path,
            "user_annotation=s" => \$user_annotation_path,
	           "file_type=s"       => \$file_type,
            );

# sanity check for input data
if (@file_query2) {
    @file_query && @file_query2 || die "Error: At least one file for each paired-end is required\n";
    @file_query == @file_query2 || die "Error: Unequal number of files for paired ends\n";
}

if (!($user_database_path)) {
    die "No reference genome was supplied\n";
}

if (@file_query < 1) {
    die "No FASTQ files were supplied\n";
}
# Sanity check for input ref. genome and annotation
unless ($user_database_path) {
  die "No reference genome was selected"
}


# Allow over-ride of system-level database path with user
if ($user_database_path) {
  $database_path = $user_database_path;
  unless (`grep \\> $database_path`) {
      die "Error: $database_path  the user supplied file is not a FASTA file";
  }
  my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
  print STDERR "STAR-indexing $name\n";
  system "mkdir index";
  my $STARp = "STAR";
  if (!($user_annotation_path)) {
    system "$STARp --runThreadN 4  --runMode genomeGenerate --genomeDir index --genomeFastaFiles $database_path";
    print "index without_annotation\n";
    $STAR_ARGS = join(" ", @ARGV);    
    }
   else{
    system "$STARp --runThreadN 4  --runMode genomeGenerate --genomeDir index --genomeFastaFiles $database_path --sjdbGTFfile $user_annotation_path";
    print "index with_annotation\n";
    $STAR_ARGS = join(" ", "--sjdbGTFfile $user_annotation_path", @ARGV);
    }
  
}

system "mkdir output; mkdir bam_output";

for my $query_file (@file_query) {
    my $second_file = shift @file_query2 if @file_query2;
    my $app  = "STAR";
    my $format = $file_type;
    chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;
	  if ($format eq 'PE') {
        my $align_command = "$app $STAR_ARGS --runThreadN 4 --genomeDir index --outReadsUnmapped Fastx --outFileNamePrefix ${basename}."." --readFilesIn $query_file $second_file --readFilesCommand gunzip -c";       
        report("Executing: $align_command\n");
        system $align_command;
        my $move_files = "mkdir output/$basename;mv Log* *STARgenome $basename*out *tab *Unmapped* output/$basename;mv *bam bam_output";
        system $move_files;
	  }
    elsif($format eq 'SE'){
	      my $align_command = "$app $STAR_ARGS --runThreadN 4 --genomeDir index --outReadsUnmapped Fastx --outFileNamePrefix ${basename}."." --readFilesIn $query_file --readFilesCommand gunzip -c ";
        report("Executing: $align_command\n");
        system $align_command;
        my $move_files = "mkdir output/$basename;mv Log* *STARgenome $basename*out *tab *Unmapped* output/$basename;mv *bam bam_output";
        system $move_files;

	}
}

sub report {
    print STDERR "$_[0]\n";
}
