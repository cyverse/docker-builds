#!/usr/bin/perl
# Downloaded from http://www.bioinformatics-made-simple.com

use strict;
use Bio::SeqIO;
my %unique;

my ($file) = @ARGV;
my $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$file.dup_removed.fasta", -format => "fasta");

while(my $seqs = $seqio->next_seq) {
  my $id  = $seqs->display_id;
  my $seq = $seqs->seq;
  unless(exists($unique{$seq})) {
    $outseq->write_seq($seqs);
    $unique{$seq} +=1;
  }
}