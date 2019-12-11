#!/usr/bin/perl -w
use strict;
use warnings;

use Getopt::Long;

my $ggFile;
my $inDir;
my $outDir;
my $debug = 1;

GetOptions("ggFile=s" => \$ggFile,
"inDir=s" => \$inDir,
"outDir=s" => \$outDir);

system "mkdir -p $outDir" unless -e $outDir;

my @taxaList;
my $taxon;
my $inLine;
my @inList;
my $element;
my %ggHoA; #key=taxonId, values = array of seqId
my %nClusterHash; #key=seqId, value=number of clusters involved (0 = unclustered)
open GG,"<$ggFile" or die "Can't open input file $ggFile: $!\n";
# read GG
while($inLine = <GG>) {
    chomp $inLine;

    @inList = split /\s+/, $inLine;

    $taxon = shift @inList;
    $taxon =~ m/(\S+):/;

    push (@taxaList, $1);

    $ggHoA{$1} = [ @inList ];
    
    foreach $element (@inList) {
        $nClusterHash{$element} = 0;
    }
    
} # read GG
close (GG);

@taxaList = sort @taxaList;

my $inIDX = $inDir .'/' . 'orthomcl.index'; #JD update
my $inMCL = $inDir . '/' . 'orthomcl.mclout'; #JD update

open IDX,"<$inIDX" or die "Can't open input file $inIDX: $!\n";
open MCL,"<$inMCL" or die "Can't open input file $inMCL: $!\n";

my %idxHash;
# key = mcl index, value = seqId 
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

my %mclHoA; # key = clusterId, values = all mcl seqIds in the cluster
my $clusterId;
my @mclSeqs;
my $mclSeq;
while($inLine = <MCL>) {
    chomp $inLine;
    
    if ($inLine =~ m/^(\d+)\s+(.*)\$/) {
        # inputline is one complete cluster (within the $ in the end)
        $clusterId = $1;
        @mclSeqs = split /\s+/, $2;
        foreach $mclSeq (@mclSeqs) {
            push @{ $mclHoA{$clusterId} }, $mclSeq;
        }
    }
    elsif ($inLine =~ m/^(\d+)\s+(.*)/) {
        # inputline is incomplete (without the $ in the end)
        $clusterId = $1;
        @mclSeqs = split /\s+/, $2;
        foreach $mclSeq (@mclSeqs) {
            push @{ $mclHoA{$clusterId} }, $mclSeq;
        }
    }
    elsif ($inLine =~ m/^\s+(.*)\$/) {
        @mclSeqs = split /\s+/, $1;
        foreach $mclSeq (@mclSeqs) {
            push @{ $mclHoA{$clusterId} }, $mclSeq;
        }
    }
    elsif ($inLine =~ m/^\s+(.*)/) {
        @mclSeqs = split /\s+/, $1;
        foreach $mclSeq (@mclSeqs) {
            push @{ $mclHoA{$clusterId} }, $mclSeq;
        }
    }
    else {
        # do nothing
    }

}# read MCL

close (MCL);

my $outIDX = $outDir . '/' . 'orthomcl.index'; #JD update
my $outMCL = $outDir . '/' . 'orthomcl.mclout'; #JD update
my $logFile = $outDir .  '/' . 'appendUnclustered.log'; #JD update

open IDX,">$outIDX" or die "Can't open output file $outIDX: $!\n";
open MCL,">$outMCL" or die "Can't open output file $outMCL: $!\n";
open LOG, ">$logFile" or die "Can't open output file $logFile: $!\n";
foreach $element (sort {$a <=> $b} keys %idxHash) {
    print IDX "$element\t$idxHash{$element}\n";
}
# for each cluster
foreach $clusterId (sort {$a <=> $b} keys %mclHoA) {
    print MCL "$clusterId\t";
    # for each seq
    foreach $element (sort {$a <=> $b} @{ $mclHoA{$clusterId} }) {
        $nClusterHash{$idxHash{$element}}++;
        print MCL "\t$element";
    }    # for each seq
    print MCL "\t\$\n";
}# for each cluster

my $nCluster = scalar (keys %mclHoA);
my $nSeq = scalar (keys %idxHash);

my %nUCHash; #key = taxon, value=nUnclustered seq
foreach $taxon (@taxaList) {
    $nUCHash{$taxon} = 0;
    foreach $element (@{ $ggHoA{$taxon} }) {
        if ($nClusterHash{$element} == 0) {
            $nUCHash{$taxon}++;
            print IDX "$nSeq\t$element\n";
            print MCL "$nCluster\t$nSeq\t\$\n";
            $nSeq++;
            $nCluster++;
        }
        elsif ($nClusterHash{$element} > 1) {
            if ($debug) {
                print "$element found in $nClusterHash{$element} clusters\n";
            }
        }
    }
}

if ($debug) {
    print LOG "inIDX = $inIDX\n";
    print LOG "inMCL = $inMCL\n";
    print LOG "\t", scalar (keys %idxHash), " seqs in ", scalar (keys %mclHoA), " clusters\n"; 
    print LOG "outIDX = $outIDX\n";
    print LOG "outMCL = $outMCL\n";
    print LOG "\t$nSeq seqs in $nCluster clusters\n"; 
    print LOG "\n";
    print LOG "Taxon\tnSeq\tnClustered\tnUnClustered\n";
    foreach $taxon (@taxaList) {
        print LOG "$taxon\t", scalar @{ $ggHoA{$taxon} }, "\t", ((scalar @{ $ggHoA{$taxon} }) - $nUCHash{$taxon}) , "\t$nUCHash{$taxon}\n";
    }
}
close IDX;
close MCL;
close LOG;

exit(0);

