#!/usr/bin/perl -w

# adapted from queryMCL.1.pl
# get the list of clusters for each taxon and taxa-pair 
# mdate 05/15/2007
#   add nSeq, nGroup to summary info (STDOUT)
# mdate 03/12/2007
#   bug fix: check if there are UniPre/UniAbs/Shared groups before output
# mdate 02/12/2007
#   add printPair opt
# mdate 11/14/2006
#   add opt
#   output cluster list for UniPre (+) and UniAbs (-) for each taxon 
# mdate 05/22/2006
#   possible to use seqId in the form of (\w+)
# cdate 11/01/2005

#Usage
#with unclustered added
#perl clusterReport.pl --ggFile=GG_Combined.txt --mclDir=unclustered_added/ --outDir=withUnclustered/
#without unclustered added
#perl clusterReport.pl --ggFile=GG_Combined.txt --mclDir=Nov_14/mcl/ --outDir=withoutUnclustered/

use strict;
use warnings;

use Getopt::Long;

my $ggFile;
my $mclDir;
my $outDir;
my $printPair = 1;
my $debug = 1;

GetOptions("ggFile=s" => \$ggFile,
"mclDir=s" => \$mclDir,
"outDir=s" => \$outDir);

system "mkdir -p $outDir" unless -e $outDir;

my $idxFile = $mclDir . '/' . 'orthomcl.index'; #JD update
my $mclFile = $mclDir . '/' . 'orthomcl.mclout'; #JD update
my $logFile = $outDir . '/' . '0_clusterReport.log'; #JD update

# declare variables
my $inLine;    # for reading input file
my @inList;
my $inItem;
my $taxon;
my @taxaList;    # hold taxon ids
my %taxonHash; # key=seqId, value=taxonId
my %ggHoA; # key=taxonId, value=array of seqId

# read GG
open GG,"<$ggFile" or die "Can't open input file $ggFile: $!\n";
while($inLine = <GG>) {
    chomp $inLine;
    @inList = split /\s+/, $inLine;
    $taxon = shift @inList;
    $taxon =~ m/(\S+):/;
    push (@taxaList, $1);
    
    foreach $inItem (@inList) {
        $taxonHash{$inItem} = $1;
        push @{ $ggHoA{$1} }, $inItem;
    }
} # read GG
close GG;

@taxaList = sort @taxaList;

# read idx
my %idxHash; # key = mcl index, value = seqId 
open IDX,"<$idxFile" or die "Can't open input file $idxFile: $!\n";
while($inLine = <IDX>) {
    $inLine =~ m/(\d+)\s+(\S+)/;
    $idxHash{$1} = $2;
}# read idx
close IDX;

# read MCL
open MCL,"<$mclFile" or die "Can't open input file $mclFile: $!\n";
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

my %orthoHoA; # HoA, key = clusterId, value = all seqId of the a cluster
my $clusterId;
my @mclSeqs;
my $mclSeq;
while($inLine = <MCL>) {
    chomp $inLine;
    
    if ($inLine =~ m/^(\d+)\s+(.*)\$/) {
        # inLine is one complete cluster (within the $ in the end)
        $clusterId = $1;
        @mclSeqs = split /\s+/, $2;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHoA{$clusterId} }, $idxHash{$mclSeq};
        }
    }
    elsif ($inLine =~ m/^(\d+)\s+(.*)/) {
        # inLine is a new cluster (incomplete, without the $ in the end)
        $clusterId = $1;
        @mclSeqs = split /\s+/, $2;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHoA{$clusterId} }, $idxHash{$mclSeq};
        }
    }
    elsif ($inLine =~ m/^\s+(.*)\$/) {
        # inLine completes a cluster
        @mclSeqs = split /\s+/, $1;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHoA{$clusterId} }, $idxHash{$mclSeq};
        }
    }
    elsif ($inLine =~ m/^\s+(.*)/) {
        # inLine continues a cluster
        @mclSeqs = split /\s+/, $1;
        foreach $mclSeq (@mclSeqs) {
            push @{ $orthoHoA{$clusterId} }, $idxHash{$mclSeq};
        }
    }
    else {
        # do nothing
    }

}# read MCL
close MCL;

my %seqCountHoH; # key1 = clusterId, key2 = taxonId, value = nSeq from the taxon in the cluster
my $seq;
my %nTaxaHash; #key = clusterId, value = nTaxa
my %nSeqHash; #key = clusterId, value = nSeq

my %uniPreHoA; # key = taxonId, value = array of clusterIds that only have seq from this taxon
my %uniAbsHoA; # key = taxonId, value = array of clusterIds that only this taxon is missing
my $taxon1;
my $taxon2;
my $pair;
my %clusterHoA; # key = taxonId1_taxonId2, value = array of clusterIds shared by the pair

# for each cluster
foreach $clusterId (sort (keys %orthoHoA)) {
    # reset seqCount
    foreach $taxon (@taxaList) {
        $seqCountHoH{$clusterId}{$taxon} = 0;
    }
    
    # for each seq
    foreach $seq (sort @{ $orthoHoA{$clusterId} }) {
        # update count
        $seqCountHoH{ $clusterId }{ $taxonHash{$seq} }++;
    }# for each seq

    # reset
    $nTaxaHash{$clusterId} = 0;
    $nSeqHash{$clusterId} = 0;
    foreach $taxon (@taxaList) {
        if ( $seqCountHoH{$clusterId}{$taxon} > 0) {
            $nTaxaHash{$clusterId}++;
            $nSeqHash{$clusterId} += $seqCountHoH{$clusterId}{$taxon};
        }    
    }

    if ($nTaxaHash{$clusterId} == 1) {
        # if only has 1 taxon => UniPre
        foreach $taxon (@taxaList) {
            if ($seqCountHoH{$clusterId}{$taxon} > 0) {
                push @{ $uniPreHoA{$taxon} }, $clusterId;
            }
        }
    }
    elsif ($nTaxaHash{$clusterId} == ((scalar @taxaList) - 1)) {
        # if has (nTaxa - 1) taxa => UniAbs
        foreach $taxon (@taxaList) {
            if ($seqCountHoH{$clusterId}{$taxon} == 0) {
                push @{ $uniAbsHoA{$taxon} }, $clusterId;
                last;
            }
        }
    }
    else {
        # do nothing
    }

    foreach $taxon1 (@taxaList) {
        foreach $taxon2 (@taxaList) {
            $pair = $taxon1 . '_' . $taxon2;
            if ( $seqCountHoH{$clusterId}{$taxon1} > 0 && $seqCountHoH{$clusterId}{$taxon2} > 0 ) {
                push @{ $clusterHoA{$pair} }, $clusterId;
            }
        }
    }
    
} # for each cluster

# output lists
my $outFile;
foreach $taxon1 (@taxaList) {
    $outFile = $outDir . '/' . $taxon1 . "+.group"; #JD update
    open OUT,">$outFile" or die "Can't open output file $outFile: $!\n";
    if (exists $uniPreHoA{$taxon1}) {
        foreach $clusterId (sort @{ $uniPreHoA{$taxon1} }) {
            print OUT "$clusterId($nTaxaHash{$clusterId}:$nSeqHash{$clusterId}";
            foreach $taxon (@taxaList) {
                print OUT ",$taxon:$seqCountHoH{$clusterId}{$taxon}"
            }
            print OUT ')';
            foreach $seq (@{ $orthoHoA{$clusterId} }) {
                print OUT "\t$seq($taxonHash{$seq})";
            }
            print OUT "\n";
        }
    }
    close OUT;
    
    $outFile = $outDir . '/' . $taxon1 . "-.group"; #JD update
    open OUT,">$outFile" or die "Can't open output file $outFile: $!\n";
    if (exists $uniAbsHoA{$taxon1}) {
        foreach $clusterId (sort @{ $uniAbsHoA{$taxon1} }) {
            print OUT "$clusterId($nTaxaHash{$clusterId}:$nSeqHash{$clusterId}";
            foreach $taxon (@taxaList) {
                print OUT ",$taxon:$seqCountHoH{$clusterId}{$taxon}"
            }
            print OUT ')';
            foreach $seq (@{ $orthoHoA{$clusterId} }) {
                print OUT "\t$seq($taxonHash{$seq})";
            }
            print OUT "\n";
        }
    }
    close OUT;
    
    if ($printPair) {
        foreach $taxon2 (@taxaList) {
            $pair = $taxon1 . '_' . $taxon2;
            $outFile = $outDir . '/' . "Shared_$pair.group"; #JD update
            open OUT,">$outFile" or die "Can't open output file $outFile: $!\n";
            if (exists $clusterHoA{$pair}) {
                foreach $clusterId (sort @{ $clusterHoA{$pair} }) {
                    print OUT "$clusterId($nTaxaHash{$clusterId}:$nSeqHash{$clusterId}";
                    foreach $taxon (@taxaList) {
                        print OUT ",$taxon:$seqCountHoH{$clusterId}{$taxon}"
                    }
                    print OUT ')';
                    foreach $seq (@{ $orthoHoA{$clusterId} }) {
                        print OUT "\t$seq($taxonHash{$seq})";
                    }
                    print OUT "\n";
                }
            }
            close OUT;
        }
    }
}

open LOG, ">$logFile" or die "Can't open output file $logFile: $!\n";
if ($debug) {
    print LOG "Input gg file contains ", scalar @taxaList, " taxa. \n";
    print LOG "Tax\tnSeq\tnGroup\t+(UPre)\t-(UAbs)\n";
    foreach $taxon (@taxaList) {
        print LOG "$taxon\t";
        print LOG scalar @{ $ggHoA{$taxon} }, "\t";
        $pair = $taxon . '_' . $taxon;
        print LOG scalar @{ $clusterHoA{$pair} }, "\t";
        if (exists $uniPreHoA{$taxon}) {
            print LOG scalar @{ $uniPreHoA{$taxon} }, "\t";
        }
        else {
            print LOG "0\t";
        }
        if (exists $uniAbsHoA{$taxon}) {
            print LOG scalar @{ $uniAbsHoA{$taxon} }, "\n";
        }
        else {
            print LOG "0\n";
        }
    }
    
    print LOG "\n";
    print LOG "Tax";
    foreach $taxon (@taxaList) {
        print LOG "\t$taxon";
    }
    print LOG "\n";
    foreach $taxon1 (@taxaList) {
        print LOG "$taxon1";
        foreach $taxon2 (@taxaList) {
            $pair = $taxon1 . '_' . $taxon2;
            if (exists $clusterHoA{$pair}) {
                print LOG "\t", scalar @{ $clusterHoA{$pair} };
            }
            else {
                print LOG "\t0";
            }
        }
        print LOG "\n";
    }
}

exit(0);

