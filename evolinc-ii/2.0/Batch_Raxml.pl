#!/usr/bin/perl
# Author: Andrew Nelson; andrew.d.l.nelson@gmail
# This is just a wrapper to infer phylogenies using RAxML in batch to keep initial script clean
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
        system("echo 'Running RAxML on $file'");
        system("raxmlHPC-PTHREADS-SSE3 --silent -s $file -n RAxML_$file -m GTRGAMMA -p 84721 -x 12345 -# 100 -f a > /dev/null 2"); #The /dev/null 2 is used to run the script silently
		system("mv RAxML_bestTree.RAxML_$file ../RAxML_families/");
		system("rm RAxML_bipartitions.RAxML_$file");
		system("rm RAxML_bipartitionsBranchLabels.RAxML_$file");
		system("rm RAxML_bootstrap.RAxML_$file");
		system("rm RAxML_info.RAxML_$file");
}
