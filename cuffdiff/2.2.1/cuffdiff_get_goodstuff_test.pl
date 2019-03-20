#!/usr/bin/perl -w
use strict;

use constant ANNOTATIONS  => '/ensembl_plant_gene_desc/';

chdir "cuffdiff_out";
system "mkdir ../sorted_data";

my (%desc);

my $fdr = shift || 0.05;

my $species = ANNOTATIONS;
open DEF, "zcat $species/*.txt.gz |" or die $!;
while (<DEF>) {
    chomp;
    my ($g,$l,$d) = split "\t";
    $g or next;
    $desc{$g}  = $d || '.';
}
close DEF;


screen_file("gene_exp.diff",1,1,0,"genes.sorted_by_fold.txt");
screen_file("gene_exp.diff",1,$fdr,0,"genes.sorted_by_fold.sig.txt");
screen_file("gene_exp.diff",0,1,0,"genes.sorted_by_expression.txt");
screen_file("gene_exp.diff",0,$fdr,0,"genes.sorted_by_expression.sig.txt");
screen_file("isoform_exp.diff",1,1,0,"transcripts.sorted_by_fold.txt",1);
screen_file("isoform_exp.diff",1,$fdr,0,"transcripts.sorted_by_fold.sig.txt",1);
screen_file("isoform_exp.diff",0,1,0,"transcripts.sorted_by_expression.txt",1);
screen_file("isoform_exp.diff",0,$fdr,0,"transcripts.sorted_by_expression.sig.txt",1);

sub screen_file {
    my $infile = shift;
    my $sortf  = shift;
    my $maxp   = shift  || 1;
    my $minfold = shift || 0;    
    my $outfile = shift;
    my $transcript = shift;

    my $index;
    if ($transcript) {
	$index = $sortf ? '-k6nr' : '-k8nr';
    }
    else {
	$index = $sortf ? '-k5nr' : '-k7nr';
    }
    
    my $header = $transcript ? "transcript\t" : '';
    $header .= join("\t",qw/gene_id gene_name sample1 sample2 fold_change
                      direction total_fpkm q-value gene_description/) . "\n";

    $outfile = "../sorted_data/$outfile";
    open OUT, ">$outfile" or die $!;
    print OUT $header;
    close OUT;
    open OUT , "| sort $index >>$outfile";

    open IN, $infile or die $!;
    my ($out,@out);
    while (<IN>) {
	next if /test_id/;
	next unless /OK/;
	my @line = split "\t";
	my $gene = $line[1];
	my $locus = $gene eq $line[2] ? '.' : $line[2];
	$gene =~ s/,\S+//;
	my $direction =$line[7] > $line[8] ? 'DOWN' : 'UP';
	my ($hi,$lo) = sort { $b <=> $a } $line[7], $line[8];
	next unless $hi && $lo;
	my $fold_change = $hi/$lo;
	next if $minfold && $fold_change < $minfold;
	my $p_val = $line[12];
	next if $maxp && $p_val > $maxp;
	my $out = $transcript ? "$line[0]\t" : '';
	$out .= join ("\t",$gene,$locus,$line[4],$line[5],sprintf("%.2f",$fold_change),$direction,sprintf("%.2f",($hi+$lo)),$p_val);
	if (defined $desc{$gene}) {
	    $out .= "\t$desc{$gene}\n";
	}
	else {
	    $out .= "\t\n";
	}

	print OUT $out;
    }
    close IN;
    close OUT;
    
}

exit 0;


