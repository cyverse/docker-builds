#!/bin/env perl

use strict;
use Pod::Usage;
use Getopt::Long;
use URI::Escape qw(%escapes);
use constant DEFAULT_SOURCE => 'vcf';
use constant DEFAULT_TYPE => 'sequence_feature';

my ($input, $output, $source, $type, $help);
GetOptions (
	"input|i=s"     => \$input,
	"output|o=s"	=> \$output,
	"source|s=s"	=> \$source,
	"type|t=s"	=> \$type,  
	"h|help"	=> \$help,
);
pod2usage(2) unless $input;
pod2usage(1) if $help;

$source||=DEFAULT_SOURCE;
$type||=DEFAULT_TYPE;

my $header_line;

if ( $input && $input ne 'stdin') {
	push @ARGV, $input;
}
if ($output) {
	open OUT, '>', $output or die "Can't open $output\n";
	select(OUT);
}

print "##gff-version 3\n";
while (<>) {
	if (/^##/) {
		print;
		if (/^##source=(.*)/) {
			if ($source eq DEFAULT_SOURCE) {
				$source=$1;
				chomp $source;
			}
		}
	} elsif (/^#CHROM/) {
		$header_line=1;
		next;
	}
	if ($header_line) {
		chomp;
		my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample)=split /\t/;
		my @attr;
		
		unless ($id eq ".") {
			push @attr, "ID=$id";
		}
		foreach my $chr (qw/% ; , = &/) {
			my $escape=$escapes{$chr};
			foreach ($filter, $info, $format, @sample) {
				s/$chr/$escape/g;
			}
		}
		push @attr, "REF=$ref", "ALT=$alt", "FILTER=$filter", "INFO=$info";
		if ($format) {
			push @attr, "FORMAT=$format";
		}
		if (@sample) {
			push @attr, "SAMPLE=" . join(",", @sample);
		}
		my @output;
		push @output, $chrom, $source, $type, $pos, $pos+length($ref)-1, $qual, ".", ".", join(";", @attr);
		print join ("\t", @output) . "\n";
	}
}

if ($output) {
	close OUT;
}

=head1 NAME

vcf_to_gff.pl - converting vcf to gff

=head1 SYNOPSIS

vcf_to_gff.pl [-i INFILE] [-o OUTFILE] [-s SOURCE] [-t TYPE]

 Options:
 [-i INPUT] = VCF input file. default is STDIN.
 [-o OUTFILE] = GFF3 output file. default is STDOUT.
 [-s SOURCE] = source column in GFF3. default is vcf or source meta field.
 [-t TYPE] = type column in GFF3. default is sequence_feature.

=head1 DESCRIPTION

This script convert vcf files to gff files.

=head1 BUGS


=head1 SUPPORT


=head1 AUTHOR

Zhenyuan Lu
Cold Spring Harbor Laboratory
luj@cshl.edu

=head1 COPYRIGHT

This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.


=head1 SEE ALSO

perl(1).

=cut

################### main pod documentation end ###################


1;
# The preceding line will help the module return a true value


1;
