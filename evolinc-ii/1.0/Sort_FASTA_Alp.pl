#!/usr/bin/perl -w

my $usage="\nUsage: $0 [-hrg] [fastaFileName1 ...]\n".
    "  -h: help\n".
    "  -r: reverse\n" .
    "  -g: remove gaps '-' from the sequence\n".
    "Sort FASTA sequences alphabetically by names.  If multiple files are \n".
    "given, sequences in all files are marged before sorting.  If no \n".
    "argument is given, it will take STDIN as the input.\n" .
    "Note that the entire sequence label including spaces is used as\n".
    "the name.\n";

our($opt_h, $opt_g, $opt_r);

use Bio::SeqIO;

use Getopt::Std;
getopts('hgr') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

my $format = "fasta";
my @seqArr = ();

@ARGV = ('-') unless @ARGV;
while (my $file = shift) {
    my $seqio_obj = Bio::SeqIO->new(-file => $file, -format => $format);
    while (my $seq = $seqio_obj->next_seq()) {
       # need to deal with spaces
	$seq->desc( $seq->id . " ". $seq->desc);

	push(@seqArr, $seq);
    }
}

if (defined($opt_r)) {
    @seqArr = sort { - ($a->desc() cmp $b->desc()) } @seqArr;
} else {
    @seqArr = sort { $a->desc() cmp $b->desc() } @seqArr;
}

my $seqOut = Bio::SeqIO->new(-fs => \*STDOUT, -format => $format);
foreach my $s (@seqArr) {
    # prints "id desc", and desc was modified, returning it to original
    my $thisDesc = $s->desc;
    $thisDesc =~ s/^\S+ //; # remove the first word.
    $s->desc($thisDesc);

    if(defined($opt_g)) {
	my $tmp = $s->seq();
	$tmp =~ s/-//g;
	$s->seq($tmp);
    }
    $seqOut->write_seq($s);
}

exit;


