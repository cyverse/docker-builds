#!/usr/bin/perl -w
my $scriptName = 'queryMCL.5.pl';
 
#Jeremy DeBarry 2/5/2015: This simple script is modified from a previous version written by the Kissinger Lab at the University of Georgia as part of a set of scripts and programs designed to cluster orthologs and assemble custom gene sets.  The source script is not licensed.  The original author is Dr. Chih Horng Kuo.  The script was modified by Jeremy DeBarry for use on iPlant Platforms.  

#usage: perl queryOrthoMCL.pl --ggFile=/path/to/gg/file --minMaxFile=/path/to/minMasFile --idxFile=/path/to/indxFile --mclFile=/path/to/mcl/file --outFile=outFileName --debug=1


use strict;
use warnings;

use Getopt::Long;
use File::Basename;

my $ggFile;
my $idxFile;
my $mclFile;
my $outFile;
my $minMaxFile; #JD ADD
my $debug;

GetOptions(    "ggFile=s"        => \$ggFile,
            "minMaxFile=s"    => \$minMaxFile, #JD ADD
            "idxFile=s"        => \$idxFile,
            "mclFile=s"        => \$mclFile,
            "outFile=s"        => \$outFile,
            "debug=i"        => \$debug);

# open input file
open GG,"<$ggFile" or die "Can't open input file $ggFile: $!\n";
open MM,"<$minMaxFile" or die "Can't open input file $minMaxFile: $!\n"; #JD add
open IDX,"<$idxFile" or die "Can't open input file $idxFile: $!\n";
open MCL,"<$mclFile" or die "Can't open input file $mclFile: $!\n";

# open output file
my $outDir = dirname($outFile);
system "mkdir -p $outDir" unless -e $outDir;
open OUT,">$outFile" or die "Can't open output file $outFile: $!\n";
my $logFile = $outDir . '/' . '0_queryOrthoMCL.log'; #JD ADD


# declare variables
my $inLine;    # for reading input file
my @inList;
my $element;
my $taxon;
my @taxaList;    # hold taxon ids
my %taxonHash; # key=seqId, value=taxonId

# read GG
while($inLine = <GG>) {
    chomp $inLine;

    @inList = split /\s+/, $inLine;

    $taxon = shift @inList;
    $taxon =~ m/(\S+):/;

    push (@taxaList, $1);

    foreach $element (@inList) {
        $taxonHash{$element} = $1;
    }
} # read GG

# close input file
close GG;

@taxaList = sort @taxaList;

# prompt users for query criteria
#my %nMaxHash; #JD REMOVE
#my %nMinHash; #JD REMOVE
#print "\n";
#print "minimum nTaxa = ";
#chomp ($nMinHash{'nTaxa'} = <STDIN>); 
#print "maximum nTaxa = ";
#chomp ($nMaxHash{'nTaxa'} = <STDIN>); 
#print "\n";
#print "minimum nSeq = ";
#chomp ($nMinHash{'nSeq'} = <STDIN>); 
#print "maximum nSeq = ";
#chomp ($nMaxHash{'nSeq'} = <STDIN>); 
#foreach $taxon (@taxaList) { #JD REMOVE
#    print "\n"; #JD REMOVE
#    print "minimum nSeq for $taxon = "; #JD REMOVE
#    chomp ($nMinHash{$taxon} = <STDIN>); #JD REMOVE
#    print "maximum nSeq for $taxon = "; #JD REMOVE
#    chomp ($nMaxHash{$taxon} = <STDIN>); #JD REMOVE
#} #JD REMOVE

#read MM #JD ADD
my %nMaxHash; #JD ADD
my %nMinHash; #JD ADD

while ($inLine = <MM>) { #JD ADD
    chomp $inLine; #JD ADD
    
    @inList = split (',', $inLine); #JD ADD
    ($nMinHash{$inList[0]} = $inList[1]); #JD ADD
    ($nMaxHash{$inList[0]} = $inList[2]); #JD ADD
} #read MM #JD ADD
close MM; #JD ADD

my %idxHash;
# key = mcl index, value = AA_SEQUENCE_ID 
# read idx
while($inLine = <IDX>) {
    $inLine =~ m/(\d+)\s+(\S+)/;
    $idxHash{$1} = $2;
}# read idx
close IDX;

# read MCL
#(mclheader
#mcltype matrix
#dimensions 14564x3644
#)
#(mclmatrix
#begin
#0        6301  6303  6310  6316  6325  6326  6329  6330  6331  6332  6338
#        11503 11508 11509 11520 11521 11523 11529 11537 11546 11547 11562
#        11563 11586 11588 11606 $
#3642    14241 14243 $
#3643    14414 14416 $
#)


# %orthoHash is a hash of arrays
#    key = AA_SEQUENCE_GROUP_ID
#    value = all AA_SEQUENCE_ID of the a group

my %orthoHash;
my $groupId;
my @mclSeqs;
my $mclSeq;
while($inLine = <MCL>) {
    chomp $inLine;
    
    if ($inLine =~ m/^(\d+)\s+(.*)\$/) {
        # inputline is one complete group (within the $ in the end)
        $groupId = $1;
        @mclSeqs = split /\s+/, $2;
        foreach $mclSeq (@mclSeqs)
        {
            push @{ $orthoHash{$groupId} }, $idxHash{$mclSeq};
        }
    }
    elsif ($inLine =~ m/^(\d+)\s+(.*)/) {
        # inputline is incomplete (without the $ in the end)
        $groupId = $1;
        @mclSeqs = split /\s+/, $2;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHash{$groupId} }, $idxHash{$mclSeq};
        }
    }
    elsif ($inLine =~ m/^\s+(.*)\$/) {
        @mclSeqs = split /\s+/, $1;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHash{$groupId} }, $idxHash{$mclSeq};
        }
    }
    elsif ($inLine =~ m/^\s+(.*)/) {
        @mclSeqs = split /\s+/, $1;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHash{$groupId} }, $idxHash{$mclSeq};
        }
    }
    else {
        # do nothing
    }
}# read MCL
close (MCL);

my $countGroup = 0;
my $countSeq = 0;

my %seqCount; # key = GGTaxonId, value = nSeq from the taxon in the group
my $seq;
my $boolPrint;
my $nTaxa;
my $nSeq;

# for each group
foreach my $groupId (sort { $a <=> $b } (keys %orthoHash)) {
    # reset $boolPrint to TRUE
    $boolPrint = 1;
    
    # reset
    $nTaxa = 0;
    $nSeq = 0;
    # reset seqCount
    foreach $taxon (@taxaList) {
        $seqCount{$taxon} = 0;
    }

    # update seqCount for each taxon
    foreach $seq (sort @{ $orthoHash{$groupId} }) {
        # update count
        $seqCount{ $taxonHash{$seq} }++;
    }# update seqCount for each taxon

    # update nTaxa and nSeq
    foreach $taxon (@taxaList) {
        if ( $seqCount{$taxon} > 0) {
            $nTaxa++;
            $nSeq += $seqCount{$taxon};
        }    
    } # update nTaxa and nSeq
    
    # check if the group satisfies query criteria
#        # 1. check nTaxa in the group
#         if ($nTaxa > $nMaxHash{'nTaxa'} || $nTaxa < $nMinHash{'nTaxa'}) {
#            # set $boolPrint to FALSE
#            $boolPrint = 0;
#         }
#         # 2. check nSeq in the group
#         elsif ($nSeq > $nMaxHash{'nSeq'} || $nSeq < $nMinHash{'nSeq'}) {
#            # set $boolPrint to FALSE
#            $boolPrint = 0;
#         }
         # 3. check nSeq for each taxon
#         else {
            foreach $taxon (@taxaList) {
                if ( $seqCount{$taxon} > $nMaxHash{$taxon} || $seqCount{$taxon} < $nMinHash{$taxon} ) {
                    # set $boolPrint to FALSE
                    $boolPrint = 0;
                    # break out of the loop, no need to check further
                    last;
                }
            }     
#         }
    # check if the group satisfies query criteria

    # if satisfied all query criteria, write to output
    if ($boolPrint == 1) {
        $countGroup++;
        # write to output
        print OUT "$groupId".'('."$nTaxa:$nSeq";
        foreach $taxon (@taxaList) {
            print OUT ",$taxon:$seqCount{$taxon}"
        }
        print OUT ')';
        
        foreach $seq (@{ $orthoHash{$groupId} }) {
            print OUT "\t$seq".'('."$taxonHash{$seq}".')';
            $countSeq++;
        }
        print OUT "\n";
    
    }     # if satisfied all query criteria, write to output

} # for each group

# close output file
close(OUT);

open LOG, ">$logFile" or die "Can't open output file $logFile: $!\n"; #JD ADD AND EACH INSTANCE OF 'LOG' BELOW
if ($debug) {
    print LOG "\n", scalar(localtime), "\n";
    print LOG "Input file (gg) = $ggFile\n";
    print LOG "Input file (mm) = $minMaxFile\n"; #JD ADD
    print LOG "Input file (idx) = $idxFile\n";
    print LOG "Input file (mcl) = $mclFile\n";
    print LOG "Output file (orthogroup) = $outFile\n";
    
    print LOG "Input gg file contains ", scalar @taxaList, " taxa. ";
    
    print LOG "Query criteria = \n";
    print LOG "\tTax\tnMin\tnMax\n";
#    print LOG "\tnTaxa\t$nMinHash{'nTaxa'}\t$nMaxHash{'nTaxa'}\n";
#    print LOG "\tnSeq\t$nMinHash{'nSeq'}\t$nMaxHash{'nSeq'}\n";
    foreach $taxon (@taxaList) {
        print LOG "\t$taxon\t$nMinHash{$taxon}\t$nMaxHash{$taxon}\n";
    }
    print LOG "\n";
    
    print LOG "The output file has $countGroup groups and $countSeq sequences.\n";
}

exit(0);

