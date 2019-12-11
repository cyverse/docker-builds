#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::SearchIO;

#Jeremy DeBarry 10/8/2014: This simple script is modified from a previous version written by the Kissinger Lab at the University of Georgia as part of a set of scripts and programs designed to cluster orthologs and assemble custom gene sets.  The source script is not licensed.  The original author is Dr. Chih Horng Kuo.  The script was modified by Jeremy DeBarry for use on iPlant Platforms.  

#Script will parse a BLAST output into a 'bpo' file format, useful as input into the OrthoMCL program.

#usage: perl parseBlastBpo.pl --inFile path/to/BLASToutput --outDir path/to/output/directory/ 

my ($inFile, $outFile, $debug);

GetOptions(	"inFile=s"	=> \$inFile,
			"outFile=s"	=> \$outFile,
			"debug=i"	=> \$debug);
			

my $searchio = new Bio::SearchIO(-format=> 'blast', -file=> "$inFile");

# open output file
open OUT,">$outFile" or die "Can't open output file $outFile\n";

#FORMAT of all.bpo or "usr_bpo_file"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1;At1g01190;535;At1g01190;535;0.0;97;1:1-535:1-535.
#2;At1g01190;535;At1g01280;510;2e-56;29;1:69-499:28-474.
#3;At1g01190;535;At1g11600;510;1e-45;27;1:59-531:21-509.
#
#Each line represents each query-subject similarity relation. And all the info is
#separated by ";", which are, in order, 
#	(1)similarity id, 
#	(2)query id, 
#	(3)query length, 
#	(4)subject id, 
#	(5)subject length, 
#	(6)BLAST E-value, 
#	(7)percent identity, 
#	(8)HSP info (each HSP is in the format of 
#		(8.1)HSP_id:
#		(8.2)query_start-query_end:
#		(8.3)subject_start-subject_end. 
#	   different HSP info are seperated by "." )
#IMPORTANT: 1. Similarity ID represents BPO file line id, so it should start 
#              from 1 for the first line, and be consecutive for the whole file.
#           2. BPO file is a parsing result from BLAST, so for each query gene
#              id, its hits can't be scattered in the file, but should be listed 
#              in ajacent lines.


my $countBpo = 0;
while( my $result = $searchio->next_result ) {
	my $queryId = $result->query_name;
	my $queryLength = $result->query_length;
	print "$queryId\t$queryLength\n" if ($debug);
	
	while( my $hit = $result->next_hit ) {
		$countBpo++;
		
		my $subjectId = $hit->name;
		my $subjectLength = $hit->length;
		print "\t$subjectId\t$subjectLength\n" if ($debug);
		my $significance = $hit->significance;
		my $matchLength = 0;
		my $nIdentical = 0;
		my $nConserved = 0;

		my $hspEntry;
		my @hspArray;

		my $countHSP = 0;
		while( my $hsp = $hit->next_hsp ) {
			$countHSP++;	
			my $nIdenticalHSP = $hsp->num_identical;
			$nIdentical = $nIdentical + $nIdenticalHSP;
			my $nConservedHSP = $hsp->num_conserved;
			$nConserved = $nConserved + $nConservedHSP;
			my $matchLengthHSP = $hsp->length('total');
			$matchLength = $matchLength + $matchLengthHSP;

			my $subjectStartHSP = $hsp->start('hit');
			my $subjectEndHSP = $hsp->end('hit');
			my $queryStartHSP = $hsp->start('query');
			my $queryEndHSP = $hsp->end('query');

			$hspEntry = "$countHSP:$queryStartHSP-$queryEndHSP:$subjectStartHSP-$subjectEndHSP.";
			push @hspArray, $hspEntry;
			
			print "\t\t$hspEntry\n" if ($debug);

		}#while hsp
		
		my $percentIdentity = 0;
		if ($matchLength > 0) {
			$percentIdentity = $nIdentical / $matchLength * 100;
		}
		print OUT "$countBpo;$queryId;$queryLength;$subjectId;$subjectLength;";
		
		printf OUT "%.2e", $significance;
		print OUT ";";
		printf OUT "%.0f", $percentIdentity;
		print OUT ";";

		foreach $hspEntry (@hspArray) {
			print OUT $hspEntry;
		}
		print OUT "\n";

	}#while subject
} #while query

close (OUT);

exit(0);


