#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper;

use constant CUFFLINKS  => ('1.3.0' => '/cufflinks-1.3.0.Linux_x86_64/',
                            '2.0.2' => '/cufflinks-2.0.2.Linux_x86_64/',
                            '2.1.1' => '/cufflinks-2.1.1.Linux_x86_64/',
                            '2.2.0' => '/cufflinks-2.2.0.Linux_x86_64/',
			    '2.2.1' => '/cufflinks-2.2.1.Linux_x86_64/');

report_input_stack();

my $cannedReference      = '';
my $userReference        = '';
my (@sampleOne,@sampleTwo,@sampleThree,@sampleFour,@sampleFive,@sampleSix,@sampleSeven,@sampleEight,@sampleNine,@sampleTen,@sampleEleven,@sampleTwelve,@sampleThirteen,@sampleFourteen,@sampleFifteen,@sampleSixteen,$fdr,$mbf);
my ($tophat_out,$cuffmerge_out,$nameOne,$nameTwo,$nameThree,$nameFour,$nameFive,$nameSix,$nameSeven,$nameEight,$nameNine,$nameTen,$nameEleven,$nameTwelve,$nameThirteen,$nameFourteen,$nameFifteen,$nameSixteen,$version,$mask_file);
my $fasta = '';
my $result = GetOptions (
                          'cannedReference=s'     => \$cannedReference,
                          'userReference=s'       => \$userReference,
                          'sampleOne=s'           => \@sampleOne,
                          'sampleTwo=s'           => \@sampleTwo,
                          'sampleThree:s'         => \@sampleThree,
                          'sampleFour:s'          => \@sampleFour,
                          'sampleFive:s'          => \@sampleFive,
                          'sampleSix:s'           => \@sampleSix,
                          'sampleSeven:s'         => \@sampleSeven,
                          'sampleEight:s'         => \@sampleEight,
                          'sampleNine:s'          => \@sampleNine,
                          'sampleTen:s'           => \@sampleTen,
                          'sampleEleven=s'        => \@sampleEleven,
                          'sampleTwelve=s'        => \@sampleTwelve,
                          'sampleThirteen:s'      => \@sampleThirteen,
                          'sampleFourteen:s'      => \@sampleFourteen,
                          'sampleFifteen:s'       => \@sampleFifteen,
                          'sampleSixteen:s'       => \@sampleSixteen,
                          'nameOne=s'             => \$nameOne,
                          'nameTwo=s'             => \$nameTwo,
                          'nameThree:s'           => \$nameThree,
                          'nameFour:s'            => \$nameFour,
                          'nameFive:s'            => \$nameFive,
                          'nameSix:s'             => \$nameSix,
                          'nameSeven:s'           => \$nameSeven,
                          'nameEight:s'           => \$nameEight,
                          'nameNine:s'            => \$nameNine,
                          'nameTen:s'             => \$nameTen,
                          'nameEleven=s'          => \$nameEleven,
                          'nameTwelve=s'          => \$nameTwelve,
                          'nameThirteen:s'        => \$nameThirteen,
                          'nameFourteen:s'        => \$nameFourteen,
                          'nameFifteen:s'         => \$nameFifteen,
                          'nameSixteen:s'         => \$nameSixteen,
			  'version=s'             => \$version,
			  'mask-file=s'           => \$mask_file,
			  'FDR=s'                 => \$fdr,
			  'frag-bias-correct=s'   => \$fasta,
                          'max-bundle-frags=s'    => \$mbf
			  
);

# Annotation sanity check
unless ($cannedReference || $userReference) {
    die "Reference or custom annotations must be supplied for CuffDiff\n";
}
# Custom trumps canned
if ($userReference) {
    $cannedReference = $userReference;
}


my %ver = CUFFLINKS;
my $cufflinksp = $ver{$version} || die "Version $version of Cufflinks is not supported\n";
my $cmd = $cufflinksp . 'cuffdiff';


# get rid of empty labels (assumes in order)
my @labels = grep {$_} ($nameOne,$nameTwo,$nameThree,$nameFour,$nameFive,$nameSix,$nameSeven,$nameEight,$nameNine,$nameTen,$nameEleven,$nameTwelve,$nameThirteen,$nameFourteen,$nameFifteen,$nameSixteen);
my $ARGS = join(' ', @ARGV);
$cmd .= " $ARGS -o cuffdiff_out";

if ($fdr) {
    $cmd .= " --FDR $fdr";
}

if ($mbf) {
    $cmd .= " --max-bundle-frags $mbf";
}

if ($fasta) {
    $cmd .= " --frag-bias-correct $fasta";
}

if ($mask_file) {
    $cmd .= " -M $mask_file";
}

if (@labels) {
    $cmd .= ' --labels '.join(',',@labels).' ';
}

$cmd .= " -p 6 ";

# Append GTF file
# but check it first for p_id and CDS features (cummeRbund will puke otherwise)
if (! `grep p_id $cannedReference` || ! `grep CDS $cannedReference`) {
    system "/cuffdiff_fix_annotations.pl $cufflinksp $cannedReference $fasta";
}

$cmd .= $cannedReference ? "$cannedReference " : ' ';


for (\@sampleOne,\@sampleTwo,\@sampleThree,\@sampleFour,\@sampleFive,\@sampleSix,\@sampleSeven,\@sampleEight,\@sampleNine,\@sampleTen,\@sampleEleven,\@sampleTwelve,\@sampleThirteen,\@sampleFourteen,\@sampleFifteen,\@sampleSixteen) {
    if ($_ && @$_ == 1) {
	if (-d $_->[0]) {
	    my $d = shift @$_;
	    while (my $f = <$d/*>) {
		push @$_, $f
		}
	}
    }
    
    $cmd .= ' '.join(',',@$_) if @$_;
}

print STDERR "Running $cmd\n";

system("$cmd 2>cuffdiff.stderr");

my $success = -e "cuffdiff_out/gene_exp.diff";

exit 1 unless $success;


# sort the data some
system "/cuffdiff_get_goodstuff_test.pl $fdr";


# Now plot a few graphs
my @pairs;
my %pair;
for my $i (@labels) {
    for my $j (@labels) {
	my $pair = join '.', sort ($i,$j);
	next if $pair{$pair}++ || $i eq $j;
	push @pairs, [$i,$j];
    }
}

my $r = <<END;
library(cummeRbund)
cuff <- readCufflinks()
png('../graphs/density_plot.png')
csDensity(genes(cuff))
dev.off()
END
;

for (@pairs) {
    my ($i,$j) = @$_;
    $r .= <<END;
png("../graphs/$i\_$j\_scatter_plot.png")
csScatter(genes(cuff),"$i","$j",smooth=T)
dev.off()
png("../graphs/$i\_$j\_volcano_plot.png")
csVolcano(genes(cuff),"$i","$j",alpha=0.05,showSignificant=T);
dev.off()
END
;

}

open RS, ">cuffdiff_out/basic_plots.R";
print RS $r;
close RS;

system "mkdir graphs";
chdir "cuffdiff_out";
system "perl -i -pe 's/(\\S)\\s+?\\#/$1\\-/' *.*";
system "R --vanilla < basic_plots.R";

system "rm -f ../munged.gtf" if -e "../munged.gtf";

exit 0;

sub report {
    print STDERR "$_[0]\n";
}


sub report_input_stack {
    my @stack = @ARGV;
    my %arg;

    while (@stack) {
        my $k = shift @stack;
    my $v = shift @stack;
        if ($v =~ /^-/) {
            unshift @stack, $v;
            $v = 'TRUE';
        }
        push @{$arg{$k}}, $v;
    }

    report("Input parameters:");
    for (sort keys %arg) {
        report(sprintf("%-25s",$_) . join(',',@{$arg{$_}}));
    }
}

