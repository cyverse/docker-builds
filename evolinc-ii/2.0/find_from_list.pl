#!/usr/bin/perl -w
#Author: Eric Lyons
use strict;
use Data::Dumper;

my ($file1, $file2) = @ARGV;

my @list;
open (IN, $file1) or die "Poops";
while (<IN>)
 {
  s/\n|\r//g;
  s/$/_/g;
  push @list, $_ if $_;
 }	

close IN;
#print Dumper \@list;


$/="\n>";
open (IN, $file2);
my %data;
while (<IN>)
 {
  s/>//g;
  foreach my $item (@list)
   {
     #print ">".$_ if /$item/;
     push @{$data{$item}},">".$_."\n" if /$item/;
   }
 }
close IN;
#print Dumper \%data;

#make files and dump goods
`mkdir -p lincRNA_families` unless -d "output";
foreach my $item (sort keys %data)
 {
   open (OUT, ">./lincRNA_families/$item.FASTA") || die "Can't open outfile: /lincRNA_families/$item: $!";
   my $string = join "\n", @{$data{$item}};
   $string =~ s/\n+/\n/g;
   print OUT $string;
   close OUT;
}

