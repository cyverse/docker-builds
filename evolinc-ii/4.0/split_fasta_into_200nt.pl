#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my ($file) = @ARGV;
my $len = 200;
my $over = 0;
my ($seq_id, $seq);
my $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");

while(my $seqs = $seqio->next_seq) {
  my $id  = $seqs->display_id;
  my $seq = $seqs->seq;

for (my $i = 1; $i <= length $seq; $i += ($len - $over)) {
    my $s = substr ($seq, $i - 1, $len);
    if (length $s > 199)
    {print "$id($i-", $i + (length $s) - 1, ")\n$s\n"}
    else {next}
}
}