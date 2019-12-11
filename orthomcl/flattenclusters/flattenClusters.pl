#!/usr/bin/perl -w
my $scriptName = 'flatOrthoGroup.2.pl';

#Jeremy DeBarry 2/5/2015: This simple script is modified from a previous version written by the Kissinger Lab at the University of Georgia as part of a set of scripts and programs designed to cluster orthologs and assemble custom gene sets.  The source script is not licensed.  The original author is Dr. Chih Horng Kuo.  The script was modified by Jeremy DeBarry for use on iPlant Platforms.  

#usage: perl flattenClusters.pl --inFile=/path/to/input --outFile=name_of_out_file --mapFile=/path/to/map/file --debug=1

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

my $inFile;
my $outFile;
my $mapFile;
my $debug;

GetOptions(    "inFile=s"        => \$inFile,
            "outFile=s"        => \$outFile,
            "mapFile=s"        => \$mapFile,
            "debug=i"        => \$debug);

my $outDir = dirname($outFile);
system "mkdir -p $outDir" unless -e $outDir;
my $logFile = $outDir . '/' . '0_flattenClusters.log'; #JD ADD

my $inLine;
my @inList;
my $inItem;

my %mapHash; #key=seqId, value=seqInfo

if ($mapFile) {
    open IN,"<$mapFile" or die "Can't open input file $mapFile: $!\n";
    while ($inLine = <IN>) {
        chomp $inLine;
        @inList = split /\t/, $inLine;
        $mapHash{$inList[0]} = $inList[1];
    }
    close IN;
}

my $groupId;
my $nSeqStat;
my $seqId;
my $taxonId;

my $countGroup = 0;
my $countSeq = 0;

open IN,"<$inFile" or die "Can't open input file $inFile: $!\n";
open OUT,">$outFile" or die "Can't open output file $outFile: $!\n";

# read input
while($inLine = <IN>) {
# input format 
# GROUP_ID(nTaxa:nSeq,ggTaxonId1:nSeqTaxon1,ggTaxonId2:nSeqTaxon2,...) SEQ_ID1(TAXON_ID1)
# e.g. 12345(2:3,ch:2,ch:1,py:0,ta:0,tp:0) 100(ch) 101(ch) 200(cp)
    $countGroup++;
    chomp $inLine;

    @inList = split /\s+/, $inLine;

    $inItem = shift @inList;
    $inItem =~ m/^(\d+)\((.+)\)/;
    $groupId = $1;
    $nSeqStat = $2;
    
    # for each seq
    foreach $inItem (@inList) {
        $countSeq++;
        $inItem =~ m/^(\S+)\((\w+)\)/;
        $seqId = $1;
        $taxonId = $2;
        
        if ($mapFile) {
            if (exists $mapHash{$seqId}) {
                print OUT "$groupId\t$taxonId\t$mapHash{$seqId}\n";
            }
            else {
                die "$seqId does not exist in mapFile\n";
            }
        }
        else {
            print OUT "$groupId\t$taxonId\t$seqId\n";
        }
    }    # for each seq
    
} 
close IN;
close OUT;

open LOG, ">$logFile" or die "Can't open output file $logFile: $!\n"; #JD ADD AND EACH INSTANCE OF 'LOG' BELOW
if ($debug) {
    print LOG "Processed $countGroup groups and $countSeq sequences.\n";
}

exit(0);

