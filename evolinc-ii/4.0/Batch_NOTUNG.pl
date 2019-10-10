#!/usr/bin/perl
# Author: Andrew Nelson; andrew.d.l.nelson@gmail
# This is just a wrapper to reconcile phylogenies using Notung in batch
use strict;
use warnings;

my $listFile = $ARGV[0];
my $species_tree = $ARGV[1];
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
        system("echo 'Running Notung on $file'");
        system("java -jar /Notung-2.8.1.7.jar -s $species_tree -g RAxML_bipartitions.RAxML_$file --root --treeoutput newick --nolosses --speciestag prefix --edgeweights name");
	system("java -jar /Notung-2.8.1.7.jar -s $species_tree -g RAxML_bipartitions.RAxML_$file.rooting.0 --rearrange --threshold 70 --treeoutput newick --speciestag prefix --log --edgeweights name --events --parsable --treestats --savepng");
}
