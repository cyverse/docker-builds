#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#Jeremy DeBarry 9/30/2014: This simple script is modified from a previous version written by the Kissinger Lab at the University of Georgia as part of a set of scripts and programs designed to cluster orthologs and assemble custom gene sets.  The source script is not licensed.  The original author is Dr. Chih Horng Kuo.  The script was modified by Jeremy DeBarry for use on iPlant Platforms.  

#usage: perl fastaRename.pl --inFile path/to/inputfasta --outDir path/to/output/directory/ --taxonID user_chosen_2letter_Identifier

my $inFile;
my $outDir;
my $taxonId;

GetOptions(	"inFile=s"		=> \$inFile,
			"outDir=s"		=> \$outDir,
			"taxonId=s"		=> \$taxonId);

system "mkdir -p $outDir" unless -e $outDir;

my $outFile = $outDir ."/". $taxonId . '.fasta';
my $mapFile = $outDir ."/". $taxonId . '.map';
my $ggFile = $outDir ."/". $taxonId . '.gg';

my $inLine;
# read inFasta
my $seqNameIn;
my $seqNameOut;
my $seq;
my %seqNameHash; #key = seqNameOut, value = seqNameIn
my %seqHash; #key = seqNameOut, value = seq
my $countSeq = 0;
{
	open IN,"<$inFile" or die "Can't open input file $inFile: $!\n";
	# redefine the record separator
	local $/ = ">";
	$inLine = <IN>; # toss the first record, which only consists of ">"
	while ($inLine = <IN>) {
		$countSeq++;
		chomp $inLine;
		($seqNameIn, $seq) = split(/\n/,$inLine,2);
		$seqNameOut = $taxonId . (sprintf "%05s", $countSeq);
		$seq =~ tr/ \t\n\r//d; # Remove whitespace
		$seqNameHash{$seqNameOut} = $seqNameIn;
		$seqHash{$seqNameOut} = $seq;
	}
	close IN;
}

open OUT,">$outFile" or die "Can't open output file $outFile\n";
open MAP,">$mapFile" or die "Can't open output file $mapFile\n";
open GG,">$ggFile" or die "Can't open output file $ggFile\n";

print GG "$taxonId:";
foreach $seqNameOut (sort keys %seqNameHash) {
	print OUT "\>$seqNameOut\n$seqHash{$seqNameOut}\n";
	print MAP "$seqNameOut\t$seqNameHash{$seqNameOut}\n";
	print GG " $seqNameOut";
}
print GG "\n";

print "Number of sequences = $countSeq \n";

close OUT;
close MAP;
close GG;

exit(0);
