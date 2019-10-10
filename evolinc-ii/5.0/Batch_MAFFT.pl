#!/usr/bin/perl
#Author: Andrew Nelson; andrew.d.l.nelson@gmail.com
#Usage: perl <Batch_MAFFT.pl> <listfile>
#Script to align sequences in batch mode
use strict;
use warnings;

my $listFile = $ARGV[0];
my @list;

open (AFILE, $listFile) or die "cannot open $listFile\n";
while (my $line = <AFILE>) {
        chomp $line;
        push @list, $line;
}
close AFILE;
#print "\n@list\n\n"; #to test the elements of the array

for (my $i=0; $i<@list; $i++) {
        my $file = $list[$i];
		my $word_count = `grep -c '>' $file`;
		if ($word_count <= 1)
		{
			system("cp $file Final_results/")
		}
		else {
        system("mafft --globalpair --maxiterate 1000 --quiet $file > Aligned_$file");
		system("mv Aligned_$file Final_results/");
		}
}
