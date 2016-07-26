#!/usr/bin/perl -w
#perl fa2oneline.pl sample.fa > out.fa
# downloaded from http://www.bioinformatics-made-simple.com
# change multiline fasta to single line fasta
# use strict;

# my $input_fasta="test_out_2/Homology_Search/Hsap.Non_identified_query_lincRNAs.fasta";
# my $output_fasta="perl_test_out.fasta";

# open(IN,"<$input_fasta") || die ("Error opening $input_fasta $!");
# open(my $OUT, ">$output_fasta") || die("Error openining $output_fasta $!");

# my $line = <IN>; 
# print $OUT $line;

# while ($line = <IN>)
# {
# chomp $line;
# if ($line=~m/^>/) { print $OUT "\n",$line,"\n" ; }
# else { print $OUT $line; }
# }

use strict;

my $input_fasta=$ARGV[0];
open(IN,"<$input_fasta") || die ("Error opening $input_fasta $!");

my $line = <IN>; 
print $line;

while ($line = <IN>)
{
chomp $line;
if ($line=~m/^>/) { print "\n",$line,"\n"; }
else { print $line; }
}

print "\n";