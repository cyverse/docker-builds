#!/usr/bin/perl -w
use strict;


my $path   = shift;
my $infile = shift or die "No infile";
my $fasta  = shift || '';
$fasta = "-s $fasta" if $fasta;

print STDERR "\n\nFixing annotation file $infile to work with cuffdiff/cummeRbund!\n\n";

unless ( `grep p_id $infile` ) {
    print STDERR "\nRunning cuffcompare\n";
    system "$path/cuffcompare -r $infile $fasta -T $infile";
    system "/munge_ids.pl <cuffcmp.combined.gtf >munged.gtf";
    system "cp munged.gtf $infile";
    system "rm -f cuffcmp* *fai";
}

unless ( `cut -f3 $infile | grep CDS` ) { 
    print STDERR "\nRunning CDS addition\n";
    system "cat $infile |sed 's/exon/CDS/' > CDS.gtf";
    system "cat $infile CDS.gtf |sort -k1,1 -k4,4n -k5,5n -k3,3r >rebuilt.gtf";
    system "cp rebuilt.gtf $infile";
    system "rm -f CDS.gtf";
}
