#!/usr/bin/perl -w

# Copyright (c) 2010
# Iowa State University
# All Rights Reserved
#
# Authors: Cheng-Ting Yeh, Sanzhen Liu, and Patrick S. Schnable
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
# following conditions are met:
#   o Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following 
#     disclaimer in the documentation and/or other materials provided with the distribution.
#   o Neither the name of the Iowa State University nor the names of its contributors may be used to endorse or promote 
#     products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CHENG-TING YEH SANZHEN LIU, or PATRICK S. SCHNABLE 
# (OR IOWA STATE UNIVERSITY) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE  POSSIBILITY OF SUCH DAMAGE.
#
# DESCRIPTION: 123SNP enables the discovery of SNPs using native alignment output files produced by bowtie, novoalign, and 
#              GSNAP. Polymorphisms detected by those programs can be tallied to discover variants found when comparing 
#              reads vs. a reference genome. This script essentially uses the results from external alignment programs 
#              and performs SNP filtering via a set of specified parameters.
#
# CITATION: 123SNP can be cited as Yeh C-T, S Liu and PS Schnable, unpublished.
#

use strict;
use warnings;
use FileHandle;
use POSIX qw(ceil floor);
use Getopt::Long;
use Time::Local;

use constant true => 1;
use constant false => 0;
use constant DEFAULT_MINIMUM_READS => 3;				# Default minimum number of non-reference reads per site
use constant DEFAULT_MINIMUM_QUALITY => 15;				# Default minimum quality
use constant DEFAULT_MAXIMUM_MISMATCHES => 3;			# Default maximum number of misamtches per read
use constant DEFAULT_MAXIMUM_ALLELES => 1;				# Default maximum of alleles to encounter per site
use constant DEFAULT_MINIMUM_ALLELE_COVERAGE => 0.80;	# Default minimum allele coverage per sample
use constant DEFAULT_MINIMUM_OVERALL_COVERAGE => 0.80;	# Default read depth coverage per SNP site (all samples)
use constant DEFAULT_MAXIMUM_SPLICE_DISTANCE => 10000;	# Maximum splice distance allowed per read (GSNAP only)
use constant DEFAULT_IGNORE => 3;						# Number of bases to ignore at the beginning and end of read
use constant DEFAULT_SOURCE => "SNP_Discovery";			# Default source to use in GFF3 file
use constant DEFAULT_FEATURE => "SNP";					# Default feature to use in GFF3 file
use constant ASCII_OFFSET => 33;						# ASCII offset for qualities in bowtie and novoalign files
use constant DEFAULT_TEMP_DIRECTORY => ".";				# Default temporary directory
use constant LIMIT => 1000;

# Generate a run ID to be used for temporary files
sub generateID {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return sprintf("%d%02d%02d.%02d%02d%02d", 1900+$year, $mon+1, $mday, $hour, $min, $sec);
} # End of sub generateID

# Trim leading and trailing white space characters
sub trim {
	my $str = $_[0];
	
	if (defined($str)) {
		$str =~ s/^(\s+|\t+)//g;
		$str =~ s/(\s+|\t+)$//g;
	} # End of if statement

	return $str;
} # End of sub trim

# Returns the minimum number from a list of integers passed as arguments

sub min {
	my @sorted = sort {$a <=> $b} @_;
	return shift(@sorted);
} # End of sub min

# Returns the maximum number from a list of integers passed as arguments

sub max {
    my @sorted = sort {$a <=> $b} @_;
    return pop(@sorted);
} # End of sub max


# Given 1 single nucleotide, returns is complement base
sub complementNT {
	my $nt = $_[0];
	$nt =~ tr/AGCT/TCGA/;
	return $nt;
} # End of sub complementNT

# Given an number as argument, formats the number by adding the comma character
# every thousands, millions, etc
sub formatNumber {
    local($_) = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;

    return $_; 
} # End of sub formatNumber

# Given total seconds signifying the module run-time, format it in hh:mm:ss format
sub formatTime {
    my $totaltime = $_[0];
    my $str;

    $str = sprintf("%02d:", $totaltime / 3600);    # total hours
    $totaltime = $totaltime % 3600;
    $str .= sprintf("%02d:", $totaltime / 60);      # total minutes
    $totaltime = $totaltime % 60; 
    $str .= sprintf("%02d", $totaltime);            # total sconds

    return $str;
} # End of sub formatTime

# Counts the number of nucleotides
sub countNucleotides {
    my ($str, @substr) = @_; 
    my $totallength = length($str);
    
    foreach my $s (@substr) {
        $str =~ s/\Q$s//g;
    } # End of for each loop

    return ($totallength - length($str));
} # End of sub countNucleotides	

# Converts fastq-sanger qualities to phred qualities
sub fastq_sanger_to_phred {
    my @qual = split(//, $_[0]);
    @qual = map { ord($_) - 33 } @qual;
    
    return @qual;
} # End of sub fastq_sanger_to_phred

# Converts phred qualities to fastq-sanger
sub phred_to_fastq_sanger {
    my @qual = @_; 
    @qual = map { chr($_ + 33) } @qual;

    return join("", @qual);
} # end of phred_to_fastq_sanger

# Evaluates bowtie polymorphisms and returns the number of mismatches and their respective
# positions relative to the reference
sub getBowtiePolymorphisms {
	my ($start, $end, $strand, $poly) = @_;
	my $length = $end - $start + 1;
	my @mismatches = split(/,/, $poly);
	my @coords;

	# Reverse mismatches order if strand eq "-"
	@mismatches = reverse(@mismatches) if ($strand eq "-");

	foreach my $m (@mismatches) {
		if ($m =~ m/^(\d+):([AGTC])>([AGTC])$/i) {
			if ($strand eq "+") {
				push(@coords, sprintf("%d:%s>%s:%s", $start + $1, $2, $1+1, $3));
			} # End of if statemnet
			else {
				push(@coords, sprintf("%d:%s>%s:%s", $end - $1, $2, $length-$1, $3));
			} # End of else statement
		} # End of if statement
	} # End of foreach statement

	return @coords;
} # End of sub getBowtiePolymorphisms

# Returns the total number of mismatches of a native format read entry

sub getNativeMismatchesCount {
	my $count = 0;

	if (defined($_[0])) {
		my @tmp = split(/\s/, $_[0]);
		$count = scalar(@tmp);
	} # End of if statemnet

	return $count;
} # End of sub getNativeMismatchesCount

# Format the progress string
sub printProgress {
	my ($old_progress, $new_progress) = @_;
	my $end_of_line = false;

	if ($new_progress =~ m/\n$/) {
		$end_of_line = true;
		chomp($new_progress);
	} # End of if statement

	if (defined($old_progress)) {
		for (my $i=0; $i < length($old_progress); $i++) {
			print STDERR sprintf("\b");	# Clear previous progress
		} # End of for loop
	} # End of if statemnet

	print STDERR sprintf("%s", $new_progress);

	# adjusting/padding text with white-space
	if (defined($old_progress) && length($new_progress) < length($old_progress)) {
		for (my $i=length($new_progress); $i < length($old_progress); $i++) {
			print STDERR sprintf(" ");
		} # End of for loop
		
		for (my $i=length($new_progress); $i < length($old_progress); $i++) {
			print STDERR sprintf("\b");
		} # End of for loop
	} # End of if statement

	if ($end_of_line) {
		print STDERR sprintf("\n");
		$new_progress .= "\n";
	} # End of if statement

	return $new_progress;
} # end of sub printProgress

# prints the header lines in temp files
sub printTempHeader {
	my $fh = $_[0];
	print $fh sprintf("# %s\n", scalar(localtime(time)));
	print $fh sprintf("#\n");
	print $fh sprintf("# Col 1 - Read Name: Unique read identifier\n");
	print $fh sprintf("# Col 2 - Program: Source alignment output program from which the read was extracted for SNP calling\n");
	print $fh sprintf("# Col 3 - Chr: Chromosome in which the read aligns uniquely\n");
	print $fh sprintf("# Col 4 - Chr Coords: Chromosomal positions where the read aligns 5' -> 3' orientation\n");
	print $fh sprintf("# Col 5 - Read Coords: Read positions where they align in reference in 5' -> 3' orientation\n");
	print $fh sprintf("# Col 6 - Strand: The original read alignment orientation from bowtie/gsnap/novoalign\n");
	print $fh sprintf("# Col 7 - Quallity: ASCII-33 encoding quality value from 5' -> 3'\n");
	print $fh sprintf("# Col 8 - N Chars: Number of N characters detected in the read sequence\n");
	print $fh sprintf("# Col 9 - Mismatches: Base position in reference and reads where the polymorphism occur\n");
	print $fh sprintf("#\n");
	print $fh sprintf("# NOTES:\n");
	print $fh sprintf("#    o Chromosome and read coordinates are adjusted to be in 5' -> 3' orientation\n");
	print $fh sprintf("#    o Quality: ASCII-33 encoding (fastq-sanger) of base qualities in 5' -> 3' orientation\n");
	print $fh sprintf("#\n");
	print $fh sprintf("# Read Name\tProgram\tChr\tChr Coords\tRead Coords\tStrand\tQuality\tN Chars\tMismatches\n");
	print $fh sprintf("#\n");
} # end of sub printTempHeader

# returns the number of insertions, deletions, splice distance, and substitutions of a alignment piece
sub getGSNAPPolymorphisms {
    my ($start, $end, $est_seq, $gen_seq, $attributes) = @_; 
    my ($ins, $del, $splice, $sub) = (0, 0, 0, 0);

    if ($attributes =~ m/ins:(\d+)/i) {
        $ins = $1; 
    } # End of if statement

    if ($attributes =~ m/del:(\d+)/i) {
        $del = $1; 
    } # end of if statement

    if ($attributes =~ m/splice_dist:(\d+)/i) {
        $splice = $1; 
    } # End of if statement

    # recount number of substitutions if GSNAP detected > 0
    if ($attributes =~ m/sub:(\d+)/i && $1 > 0) {
        # compute total number of substitutions
        my @est_nucleotides = split(//, uc($est_seq));
        my @gen_nucleotides = split(//, uc($gen_seq));
        for (my $i=$start; $i <= $end; $i++) {
            if ($est_nucleotides[$i - 1] ne $gen_nucleotides[$i - 1] &&

                $est_nucleotides[$i - 1] ne "N" &&
                $gen_nucleotides[$i - 1] ne "N") {
                $sub++;
            } # End of if statement
        } # End of for loop
    } # End of if statement

    return ($ins, $del, $splice, $sub);
} # End of sub getGSNAPPolymorphisms

# returns true if splice translocation, splice inversion false otherwise

sub GSNAPSpliceTranslocationInversionScramble {
	my $translocation_inversion_scramble = false;

	if ($_[0] =~ m/splice_(translocation|inversion|scramble)/i) {

		$translocation_inversion_scramble = true;
	} # End of if statement

	return $translocation_inversion_scramble;
} # End of sub GSNAPSpliceTranslocationInversionScramble

# returns the EST coordinates of the alignment
sub getGSNAPESTCoords {
    my $coords = $_[0];
    my ($start, $end);

    if ($coords =~ m/^(\d+)\.\.(\d+)$/) {
        $start = $1;
		$end = $2;
    } # End of if statement
    else {
        print STDERR sprintf("\n\n");
        print STDERR sprintf("ERROR: INVALID EST COORDINATES: %s\n", $coords);
        print STDERR sprintf("\n");
        exit();
    } # End of else statement

    return ($start, $end);
} # End of sub getGSNAPESTCoords

# returns the Genomic coordinates of the alignment
sub getGSNAPGenomicCoords {
    my $coords = $_[0];
    my ($strand, $genomic, $start, $end);

    if ($coords =~ m/^(\+|-)(\w+):(\d+)\.\.(\d+)$/) {
        $strand = $1;
        $genomic = $2;
        $start = $3;
		$end = $4;
    } # End of if statement
    else {
        print STDERR sprintf("\n\n");
        print STDERR sprintf("ERROR: INVALID GENOMIC COORDINATES: %s\n", $coords);
        print STDERR sprintf("\n");
        exit();
    } # End of else statement

    return ($strand, $genomic, $start, $end);
} # End of sub getGSNAPGenomicCoords

# Evaluates the accumulated GSNAP paths based on the total allowed mismatches and returns
# the path number if unique and valid
sub evaluateGSNAPPaths {
	my ($paths, $paths_num, $mis_allowed, $maxSplice) = @_;
	my ($good_path, %tally);	# Keep track of the number of mismatches per path

    # filter valid reads having ins+del+subs+tail_front+tail_end <= mis_allowed
    for (my $i=1; $i <= $paths_num; $i++) {
        my $total_mismatches = $paths->{$i}->{"INSERTIONS"} +
                               $paths->{$i}->{"DELETIONS"} +
                               $paths->{$i}->{"SUBSTITUTIONS"} +
                               $paths->{$i}->{"EST_TAIL_FRONT"} +
                               $paths->{$i}->{"EST_TAIL_END"};
        
		# Check if all pieces are ok and read contains <= mis_allowed mismatches. If read
		# constains either a splice_translocation, splice_inversion, and/or splice_scramble,
		# consider the read/path to be invalid
		if ($total_mismatches <= $mis_allowed && $paths->{$i}->{"MAXIMUM_SPLICE_DISTANCE"} <= $maxSplice &&
		    $paths->{$i}->{"TRANSLOCATION_INVERSION_SCRAMBLE"} == false) {
            if (exists $tally{$total_mismatches}) {
                $tally{$total_mismatches} .= sprintf(",%d", $i);
            } # End of if statemnet
            else {
                $tally{$total_mismatches} = $i; 
            } # End of else statement
        } # End of if statement
    } # End of for loop

    # Tally contains the number of mismathes per read, see if there is a unique alignment of this read
    my @keys = sort {$a <=> $b} keys %tally;
	if (scalar(@keys) == 1) {
		my @tokens = split(/,/, $tally{$keys[0]});
		if (scalar(@tokens) == 1) {	# Ok to go
			$good_path = shift(@tokens);
		} # End of if statement
	} # End of if statement

	return $good_path;
} # End of sub evaluateGSNAPPaths

# computes mismatches positions relative to reference and read and fixes the read coordinates
# to be always facing 5' -> 3'
sub fixGSNAPPath {
	my ($name, $seq, $good_path, $paths) = @_;
	my @mismatches;

	my @est_nucleotides = split(//, $seq);
	my @est_coords = split(/\s/, $paths->{$good_path}->{"EST_COORDS"});
	my @gen_coords = split(/\s/, $paths->{$good_path}->{"GENOMIC_COORDS"});
	my @gen_sequences = split(/\s/, $paths->{$good_path}->{"GENOMIC_SEQUENCES"});
	my $strand = $paths->{$good_path}->{"STRAND"};

	# convert read sequence and all gen_sequences to uppercase
	$seq = uc($seq);
	@gen_sequences = map { uc($_) } @gen_sequences;

	for (my $i=0; $i < scalar(@est_coords); $i++) {
		my ($est_start, $est_end) = split(/\.\./, $est_coords[$i]);
		my ($gen_start, $gen_end) = split(/\.\./, $gen_coords[$i]);
		my @gen_nucleotides = split(//, $gen_sequences[$i]);


		my ($gen_pos, $gen_nt, $read_pos, $read_nt);
		for (my $offset=0; $offset < $est_end - $est_start + 1; $offset++) {
			my $read_nt = $est_nucleotides[$est_start + $offset - 1];
			my $gen_nt = $gen_nucleotides[$est_start + $offset - 1];

			# nucleotides are already in uppercase 
			if ($read_nt ne $gen_nt && $read_nt ne "N" && $gen_nt ne "N") {
				if ($strand eq "-") {	# need to fix 5' -> 3'
					$read_pos = length($seq) - $est_start - $offset + 1;
					$gen_pos = $gen_start - $offset;
					$read_nt = &complementNT($read_nt);
					$gen_nt = &complementNT($gen_nt);
				} # End of if statemnet
				else {
					$read_pos = $est_start + $offset;
					$gen_pos = $gen_start + $offset;
				} # End of if statement

				push(@mismatches, sprintf("%s:%s>%s:%s", $gen_pos, $gen_nt, $read_pos, $read_nt));
			} # End of if statement
		} # End of for loop

		if ($strand eq "-") {
			$est_coords[$i] = sprintf("%d..%d", length($seq) - $est_end + 1, length($seq) - $est_start + 1);
			$gen_coords[$i] = sprintf("%d..%d", $gen_end, $gen_start);
		} # End of if statement
	} # End of for loop

	if ($strand eq "-") {	# Keep coordinates ascending relative to reference
		@gen_coords = reverse(@gen_coords);
		@est_coords = reverse(@est_coords);
		@mismatches = reverse(@mismatches);
	} # end of if statement

	return (join(" ", @gen_coords), join(" ", @est_coords), $strand, join(" ", @mismatches));
} # End of sub fixGSNAPPath

# Given a file handler of quality file as argument, read and find the next available quality,
# returns the file handler, quality name, and quality values
sub nextGSNAPQuality {
    my $fh = $_[0];     # File handler pointer
    my ($name, $qual);
    my $stop = false;
    my $line;

    my $pos = tell($fh);        # Obtain current file position
    while (!$stop && !eof($fh)) {
        $line = <$fh>;  # Read a line
        chomp($line);   # Remove end of line character
        $line = &trim($line);   # Remove leading and trailing white space characters

        if (length($line) != 0) {
            if ($line =~ m/^>(\S+)/) {
                if (defined($name)) {
                    seek($fh, $pos, 0); # Restore position to the beginning of this line
                    $stop = true;
                } # End of if statement
                else {
                    $name = $1;
                    $qual = "";
                } # End of else statemenet
            } # End of if statement
            else {
                # Accumualte quality values
                if (length($qual) == 0) {
                    $qual = $line;
                } # end of if statement
                else {
                    $qual .= sprintf(" %s", $line);
                } # End of else statement
            } # End of else statement
        } # End of if statement

        $pos = tell($fh);
    } # End of while loop

    if (defined($name) && defined($qual)) { # Perform quality cleanup
        $qual = &trim($qual);       # Remove any leading and trailing white space characters
        my $index = index($qual, "  ");     # Find the first instance with 2 spaces
        while ($index >= 0) {
            $qual =~ s/  / /g;   # Replace all 2 spaces by 1 space
            $index = index($qual, "  ");
        } # End of while loop
    } # End of if statement
	
	if (length($qual) == 0) {	# Empty qualities
		print STDERR sprintf("\n\n");
		print STDERR sprintf("ERROR: No quality values detected for '%s'\n", $name);
		print STDERR sprintf("\n");
		exit();
	} # End of if statement

    return ($name, split(/\s/, $qual));
} # End of sub nextGSNAPQuality

# Given a read ID, quality pointer for random access and quality file handler, scan thru the entire file
# to obtain the read quality, for each quality scanned not equal to $name, update the file position in qualPTR
# to be used later
sub fetchGSNAPQuality {
	my ($name, $qualPTR, $qualFH) = @_;
	my ($filepos, $read_name, @quality);
	my $stop = false;

	if (exists $qualPTR->{$name}) {		# Quality already was read before
		# Obtain current file position before reading quality
		$filepos = tell($qualFH);

		# Since the quality we are looking for was already read, adjust the file pointer to the right
		# position and read quality
		seek($qualFH, $qualPTR->{$name}, 0);					# Adjusting
		($read_name, @quality) = &nextGSNAPQuality($qualFH);	# Read quality
		
		# Restore the file position to where its supposed to be before reading quality
		seek($qualFH, $filepos, 0);
		
		# Delete pointer since reads are supposed to be uniquely mapped, we will not need it anymore
		delete $qualPTR->{$name};
	} # End of if statement
	else {		# seek and store file positions
		while (!$stop && !eof($qualFH)) {
			$filepos = tell($qualFH);	# position in file before reading quality
			($read_name, @quality) = &nextGSNAPQuality($qualFH);

			if ($read_name eq $name) {	# Found it
				$stop = true;
			} # End of if statement
			else {	# save file position and read next quality the next iteration
				$qualPTR->{$read_name} = $filepos;
			} # End of else statement
		} # End of while loop
	} # End of else statement

	return ($qualPTR, @quality);
} # End of sub fetchGSNAPQuality

# Returns the total number of insertions and deletions before a given index in the mismatches array
sub getNovoalignInsertionsDeletionsBeforeIndex {
	my ($index, @mismatches) = @_;
	my $deletions = 0;
	my $insertions = 0;

	for (my $i=0; $i < $index; $i++) {
		if ($mismatches[$i] =~ m/^(\d+)-(\S+)$/i) {
			$deletions += length($2);
		} # End of if statement
		elsif ($mismatches[$i] =~ m/^(\d+)\+(\S+)$/i) {
			$insertions += length($2);
		} # end of else if statement
	} # end of for loop
	
	return ($insertions, $deletions);
} # End of sub getNovoalignInsertionsDeletionsBeforeIndex

# returns the number of mismatches from
sub countNovoalignMismatches {
	my $mis = $_[0];
	my $mismatches = 0;

	if (!defined($mis)) {
		$mismatches = 0;
	} # End of if statement
	else {
    	my @polymorphisms = split(/\s/, $mis);

	    foreach my $p (@polymorphisms) {
    	    if ($p =~ m/^(\d+)(.)>(.)$/i) { # mismatch polymorphism
        	    $mismatches++;   # Increment by 1
        	} # End of if statemnet
        	elsif ($p =~ m/^(\d+)\+(\S+)$/i) {  # Insertion
            	$mismatches += length($2);
        	} # End of elsif statement
        	elsif ($p =~ m/^(\d+)-(\S+)$/) {
            	$mismatches += length($2);
        	} # End of else statement
    	} # End of for each statement
	} # End of else statement

    return $mismatches;
} # End of sub countNovoalignMismatches

# Returns an array of polymorphisms computed from novoalign mismatches column. The polymorphism
# coords are from 5' -> 3' and insertions/deletions have been adjusted accordingly
sub getNovoalignPolymorphismsCoords {
	my ($ref_start, $mis) = @_;
	my @coords;

	if (defined($mis)) {
		my @mismatches = split(/\s/, $mis);
		for (my $i=0; $i < scalar(@mismatches); $i++) {
			my $poly = $mismatches[$i];
			
			if ($poly =~ m/^(\d+)(.)>(.)$/i) {  # mismatch polymorphism
				my $read_pos = $1;
				my $ref_allele = uc($2);
				my $read_allele = uc($3);
				my $ref_pos = $ref_start + $read_pos - 1;       # Position in reference of which polymorphism occurred

				# Adjusting read position based on insertions/deletions and alignment orientation
				# 9G>C 22C>G 70-C 76C>G
				# for position 76, the position is actually 75 because there is a read deletion at 70
				my ($read_insertions, $read_deletions) = &getNovoalignInsertionsDeletionsBeforeIndex($i, @mismatches);
				my $adj_read_pos = $read_pos + $read_insertions - $read_deletions;

				push(@coords, sprintf("%s:%s>%s:%s", $ref_pos, $ref_allele, $adj_read_pos, $read_allele));
			} # End of if statement
		} # End of for loop
	} # End of if statement

	return @coords;
} # End of sub getNovoalignPolymorphismsCoords

# Adjust reference and read coordinates for the regions that the read actually aligns
# if there are insertions/deletions
sub adjustNovoalignCoordinates {
	my ($ref_start, $ref_end, $read_start, $read_end, $mis) = @_;

	my (@ref_coords, @read_coords);

	# No insertion/deletion polymorphisms, leave coordinates unchanged
	if (!defined($mis) || ($mis !~ m/(\d+-\S+)/ && $mis !~ m/(\d+\+\S+)/)) {
		push(@ref_coords, sprintf("%d..%d", $ref_start, $ref_end));
		push(@read_coords, sprintf("%s..%s", $read_start, $read_end));
	} # end of if statement
	else {
		my ($ref_piece_start, $ref_piece_end) = ($ref_start, $ref_end);
		my ($read_piece_start, $read_piece_end) = ($read_start, $read_end);
		my @poly = split(/\s/, $mis);
		for (my $i=0; $i < scalar(@poly); $i++) {
			if ($poly[$i] =~ m/(\d+)(-)(\S+)/ || $poly[$i] =~ m/(\d+)(\+)(\S+)/) {	# Insertion/deletions in read
				my ($read_insertions, $read_deletions) = &getNovoalignInsertionsDeletionsBeforeIndex($i, @poly);
				my $adj_read_pos = $1 + $read_insertions - $read_deletions - 1;	# Actual read position before polymorphism
				my $length = $adj_read_pos - $read_piece_start;
				
				push(@ref_coords, sprintf("%d..%d", $ref_piece_start, $ref_piece_start + $length));
				push(@read_coords, sprintf("%d..%d", $read_piece_start, $read_piece_start + $length));
				
				if ($2 eq "-") {	# read contains missing bases (deletions)
					$ref_piece_start += $length + length($3) + 1;
					$ref_piece_end += length($3);

					$read_piece_start += $length + 1;
				} # End of if statement
				else {	# read has extra bases (insertions), expand reference
					$ref_piece_start += $length + 1;
					$ref_piece_end -= length($3);

					$read_piece_start += $length + length($3) + 1;
				} # End of else statemnet
			} # End of if statement
		} # End of for each statement
		
		push(@ref_coords, sprintf("%d..%d", $ref_piece_start, $ref_piece_end));
		push(@read_coords, sprintf("%d..%d", $read_piece_start, $read_piece_end));
	} # End of if statement

	return (join(" ", @ref_coords), join(" ", @read_coords));
} # End of adjustNovoalignCoordinates

# Returns the start and end positions of the read in the chromosome in addition to the total number of mismatches
# for read stacking filtering
sub getReadAttributes {
	my ($ref_coords, $mis) = @_;
	my ($start, $end, $poly);

	# Coordinates calculation
	my @tmp = split(/[\.\s]/, $ref_coords);
	$start = shift(@tmp);
	$end = pop(@tmp);

	# mismatches calculation
	if (defined($mis)) {
		my @tmp = split(/\s/, $mis);
		$poly = scalar(@tmp);
	} # End of if statemnet
	else {
		$poly = 0;
	} # End of else statement

	return ($start, $end, $poly);
} # end of sub getReadAttributes

# Returns the adjusted reference coordinates for trimmed reads by adding back the trimmed nucleotides 
sub getReadAttributesBeforeTrimming {
	my ($name, $ref_coords, $strand, $mis, $tlPTR, $logHandle) = @_;
	my ($filepos, $line, $read_name, $origLength, $clip_left, $clip_right, $trimmedLength, $left, $right);
	my ($ref_start, $ref_end, $mismatches);
	my $stop = false;

	if (exists $tlPTR->{$name}) {	# Trimmed coordinates already read, load it and compute
		# Obtain current file position before reading line
		$filepos = tell($logHandle);

		# Adjust file handler to the line where it will beging to read
		seek($logHandle, $tlPTR->{$name}, 0);
		$line = <$logHandle>;	# Read current line
		chomp($line);		# remove end of line character


		my ($read_name, $origLength, $clip_left, $clip_right, $trimmedLength) = split(/\t/, $line);
		$left = $clip_left - 1;
		$right = $origLength - $clip_right;		# left and right are actuall bases removed during trimming
	
		# restore file position
		seek($logHandle, $filepos, 0);
		
		# delete pointer, no need it anymore
		delete $tlPTR->{$name};
	} # End of if statmenet
	else {
		$filepos = tell($logHandle);
		while (!$stop && !eof($logHandle)) {
			$filepos = tell($logHandle);
			$line = <$logHandle>;	# Read line
			chomp($line);

			if (length($line) != 0 && $line !~ m/^#/) {
				my ($read_name, $origLength, $clip_left, $clip_right, $trimmedLength) = split(/\t/, $line);

				if ($read_name eq $name) { # Found what we needed
					$left = $clip_left - 1;
					$right = $origLength - $clip_right;
					$stop = true;
				} # End of if statemnet
				else {	# Save
					$tlPTR->{$read_name} = $filepos;
				} # End of else statement
			} # End of if statemnet
		} # End of while loop
	} # End of else statement

	if (defined($left) && defined($right)) {	# Found the bases that need to be added at the beggining and end
		($ref_start, $ref_end, $mismatches) = &getReadAttributes($ref_coords, $mis);
		
		# Adjust coordinates
		if ($strand eq "-") {
			$ref_start -= $right;
			$ref_end += $left;
		} # End of if statement
		else {
			$ref_start -= $left;
			$ref_end += $right;
		} # End of else statement
	} # End of else statemnet
	else {
		print STDERR sprintf("\n\n");
		print STDERR sprintf("ERROR: Cannot find trimmed coordinates for read '$name'\n");
		print STDERR sprintf("\n");
		exit();
	} # End of else statement

	return ($ref_start, $ref_end, $mismatches, $tlPTR);
} # end of sub getReadAttributesBeforeTrimming

# returns an array containing possible read alleles
sub getPossibleReadAlleles {
	my @tmp = split(//, $_[0]);
	my %alleles;
	foreach my $t (@tmp) {
		$alleles{$t} = 1;
	} # end of for each statement

	return keys %alleles;
} # End of sub getPossibleReadAlleles

# Given a string containing all possible alleles per site, rank the alleles based on their occurrence
sub rankAlleles {
	my ($alleles, $total_alleles) = @_;
	my @possible_alleles = &getPossibleReadAlleles($alleles);
	my %frequencies;

	# Tally based on allele frequency
	foreach my $pa (@possible_alleles) {
		my $count = &countNucleotides($alleles, $pa);
		if (exists $frequencies{$count}) {
			$frequencies{$count} .= sprintf(",%s", $pa);

		} # End of if statement
		else {
			$frequencies{$count} = $pa;
		} # End of else statement
	} # End of for each statement

	# construct sorted alleles array in descending order and return

	my @rank;
	foreach my $f (sort {$b <=> $a} keys %frequencies) {
		my @read_alleles = split(/,/, $frequencies{$f});
		my $value;
		foreach my $ra (@read_alleles) {
			for (my $i=1; $i <= $f; $i++) {
				if ($i==1) {
					$value = $ra;
				} # End of if statemetn
				else {
					$value .= $ra;
				} # End of else statemnet
			} # End of for loop
		
			push(@rank, $value);
		} # End of for loop
	} # End of for each statement

	return splice(@rank, 0, $total_alleles);
} # End of rankAlleles

# Given the reference allele and all possible allele groups, determine how many non-ref alleles are present
sub countNonReferenceAlelles {
	my ($ref_allele, @allele_groups) = @_;
	my $count = 0;

	foreach my $ag (@allele_groups) {
		$count++ if ($ref_allele ne substr($ag, 0, 1));
	} # End of for each statement

	return $count;
} # End of sub countNonReferenceAlleles

# ------------------------- MAIN PROGRAM STARTS HERE ----------------------- #
# processing command line arguments
my (@bowtie, @gsnap, $gsnapQuality, @novoalign, @native, $output);
my ($total_alleles, $num, $allele_coverage, $overall_coverage, @mismatches, $ignore, $min_quality, $assignQual, $maxSplice);
my ($source, $feature, $stacking, $tlfile, $clean, $temp);

# %unique_reads variable is to keep track of read names present in bowtie, gsnap, and novoalign output files
# in case a read was aligned more than once using different alignment programs. This variable is crucial
# for not double counting the total number of reads present per snp site

my (%unique_reads);

my $result = &GetOptions("bowtie:s{0,}" => \@bowtie,
                         "gsnap:s{0,}" => \@gsnap,
						 "gsnapQuality|gsnapQual:s{1}" => \$gsnapQuality,
						 "novoalign:s{0,}" => \@novoalign,
						 "native:s{0,}" => \@native,
						 "output|o=s{1}" => \$output,
						 "alleles|a:i{1}" => \$total_alleles,
						 "mismatches|mm|m:i{1,2}" => \@mismatches,
						 "number|n:i{1}" => \$num,
						 "allelecoverage|ac:f{1}" => \$allele_coverage,
						 "converage|c:f{1}" => \$overall_coverage,
						 "quality|q:i{1}" => \$min_quality,
						 "assignQual:i{1}" => \$assignQual,
						 "ignore:i{1}" => \$ignore,
						 "source|s:s{1}" => \$source,
						 "feature|f:s{1}" => \$feature,
						 "stacking!" => \$stacking,
						 "maxsplice:i{1}" => \$maxSplice,
						 "trimmedlog|trlog:s{1}" => \$tlfile,
						 "temp:s{1}" => \$temp,
						 "clean!" => \$clean);


unless ($result && scalar(@bowtie) + scalar(@gsnap) + scalar(@novoalign) + scalar(@native) > 0 && defined($output)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("SNP DISCOVERY USING BOWTIE, GSNAP, AND NOVOALIGN OUTPUT FILES\n");
	print STDERR sprintf("=============================================================\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("  perl %s [--bowtie <bowtie files>] [--gsnap <gsnap files>] [--gsnapQuality <gsnap quality file>]\n", $0);

	print STDERR sprintf("          [--novoalign <novoalign files>] [--native <native files] --output|-o <output.gff3> [OPTIONS]\n", $0);
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("  --bowtie <bowtie files>               : Path to alignment output files generated by bowtie. The specified\n");

    print STDERR sprintf("                                          files are assumed to only include the best and unique alignments.\n");

	print STDERR sprintf("                                          Refer to suggested bowtie program parameters below.\n");
	print STDERR sprintf("  --gsnap <gsnap files>                 : Path to alignment output files generated by gsnap. Refer to suggested\n");
	print STDERR sprintf("                                          gsnap program parameters below.\n");
	print STDERR sprintf("  --gsnapQual <gsnap quality file>      : Path to quality file containing nucleotide values.\n");
	print STDERR sprintf("                                          This option is only required if --gsnap is specified\n");
	print STDERR sprintf("                                          and --assignQual is omitted.\n");
	print STDERR sprintf("  --novoalign <novoalign files>         : Path to alignment output files generated by novoalign. Refer to\n");
	print STDERR sprintf("                                          suggested novoalign parameters below.\n");
	print STDERR sprintf("  --native <native files>               : Path to previously generated native files.\n");
	print STDERR sprintf("  --output|-o <output.gff3>             : Path to output file where SNP calls will be saved\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("  --alleles|-a <num>                    : Maximum number of read alleles to encounter per site [DEFAULT: %d]\n", DEFAULT_MAXIMUM_ALLELES);
	print STDERR sprintf("  --mismatches|-mm|-m <num>             : Maximum number of mismatches allowed per read. This option\n");
	print STDERR sprintf("                                          evaluates each read and accepts/discards based on the number\n");
	print STDERR sprintf("                                          polymorphisms regardless of read length [DEFAULT: %d]\n", DEFAULT_MAXIMUM_MISMATCHES);
	print STDERR sprintf("  --mismatches|-mm|-m <num> <length>    : Maximum number of mismatches per nucleotides of the read. This\n");
	print STDERR sprintf("                                          option is used if the input reads have variable lengths\n");
	print STDERR sprintf("                                          to allow variable tolerance in evaluating reads to be accepted\n");
	print STDERR sprintf("                                          or discared. For instance, if '--mismatches 2 36' is specified\n");
	print STDERR sprintf("                                          and the read is 75 bp in length, the total allowed number of\n");
	print STDERR sprintf("                                          mismatches of this read is: CEILING((75 x 2) / 36) = 5 mismatches\n");

	print STDERR sprintf("  --num|-n <number of reads>            : Minimum number of reads per allele in each SNP site [DEFAULT: %d]\n", DEFAULT_MINIMUM_READS);
	print STDERR sprintf("  --allelecoverage|-ac <coverage>       : Minimum allele coverage allowed for a SNP [DEFAULT: %2.2f]\n", DEFAULT_MINIMUM_ALLELE_COVERAGE);
	print STDERR sprintf("  --coverage|-c <coverage>              : Minimum overall coverage of all allowed alleles\n");
	print STDERR sprintf("                                          per SNP site [DEFAULT: %2.2f]\n", DEFAULT_MINIMUM_OVERALL_COVERAGE);
	print STDERR sprintf("  --quality|-q <quality value>          : Minimum quality value of nucleotides to consider in\n");
	print STDERR sprintf("                                          phred scale (0 ~ 40) [DEFAULT: %d]\n", DEFAULT_MINIMUM_QUALITY);
	print STDERR sprintf("  --assignQual <quality value>          : Overwrite/Assign to all nucleotides aligned by bowtie, gsnap,\n");
	print STDERR sprintf("                                          and/or novoalign by the specified quality value. Specified quality\n");

	print STDERR sprintf("                                          value is an interger value in phred scale (0 ~ 40)\n");
	print STDERR sprintf("  --ignore <number of bases>            : Specify the number of bases to ignore at the beginning and\n");
	print STDERR sprintf("                                          end of the read. Ignored nucleotides does not participate\n");
	print STDERR sprintf("                                          nor affect SNP discovery procedures [DEFAULT: %d]\n", DEFAULT_IGNORE);

	print STDERR sprintf("  --maxsplice <max. splice distance>    : Specify the maximum number of bases to allow as splice distance\n");
	print STDERR sprintf("                                          for each read (GSNAP output reads only) [DEFAULT: %s]\n", DEFAULT_MAXIMUM_SPLICE_DISTANCE);
	print STDERR sprintf("  --source <string>                     : Source string to be included in GFF3 file [DEFAULT: %s]\n", DEFAULT_SOURCE);
	print STDERR sprintf("  --feature <string>                    : Feature string to be included in GFF3 file [DEFAULT: %s]\n", DEFAULT_FEATURE);
	print STDERR sprintf("  --stacking|--nostacking               : Filters read stacking by comparing start/end chromosomal coordinates of\n");
	print STDERR sprintf("                                          of each read to determine reads stacking in certain regions for\n");
	print STDERR sprintf("                                          removal. [DEFAULT: --nostacking]\n");
	print STDERR sprintf("  --trimmedlog|trlog <trim. log file>   : Specify the trimmed log file produced by in-house fastq trimming\n");
	print STDERR sprintf("                                          script. The presence of this file will allow the program to determine\n");
	print STDERR sprintf("                                          the exact coordinates of each read when performing reads stacking removal.\n");
	print STDERR sprintf("  --temp <temp directory>               : Path to a temporary directory to save files to. Default is to save\n");

	print STDERR sprintf("                                          temporary files in current directory.\n");
	print STDERR sprintf("  --clean|--noclean                     : Enable/Disable the cleaning of temporary files [DEFAULT: --clean]\n");

	print STDERR sprintf("\n");
	print STDERR sprintf("SUGGESTED BOWTIE, GSNAP, AND NOVOALIGN PARAMETERS:\n");
	print STDERR sprintf("  %% bowtie <ebwt> -q <input.fq> --solexa-quals -a -n 2 -l 25 -m 1 -B 1 --best --strata\n");
	print STDERR sprintf("  %% gsnap -D <db_dir> -d <db_name> -B 2 -m 10 -i 2 -N 1 -n 3 <input.fas or input.fq>\n");
	print STDERR sprintf("  %% novoalign -d <index db> -f <input.fq> -F ILMFQ -R 0 -r None\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("NOTES: Change to the appropriate fastq-variant in bowtie and novoalign based on your input.\n");
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

# assigning default values
if (scalar(@mismatches) == 0) {
	push(@mismatches, DEFAULT_MAXIMUM_MISMATCHES);
} # End of if statement

$total_alleles = DEFAULT_MAXIMUM_ALLELES if (!defined($total_alleles) || $total_alleles !~ m/^\d+$/);
$num = DEFAULT_MINIMUM_READS if (!defined($num) || $num !~ m/^\d+$/);

$allele_coverage = DEFAULT_MINIMUM_ALLELE_COVERAGE if (!defined($allele_coverage) || $allele_coverage !~ m/^\d+\.?\d{0,}$/);

$overall_coverage = DEFAULT_MINIMUM_OVERALL_COVERAGE if (!defined($overall_coverage) || $overall_coverage !~ m/^\d+\.?\d{0,}$/);

$min_quality = DEFAULT_MINIMUM_QUALITY if (!defined($min_quality) || $min_quality !~ m/^\d+$/);
$ignore = DEFAULT_IGNORE if (!defined($ignore) || $ignore !~ m/^\d+$/);
$source = DEFAULT_SOURCE if (!defined($source));
$feature = DEFAULT_FEATURE if (!defined($feature));
$stacking = false if (!defined($stacking));
$maxSplice = DEFAULT_MAXIMUM_SPLICE_DISTANCE if (!defined($maxSplice) || $maxSplice !~ m/^(\d+)$/);
$clean = true if (!defined($clean));
$temp = DEFAULT_TEMP_DIRECTORY if (!defined($temp));

my $total_input = scalar(@bowtie) + scalar(@gsnap) + scalar(@novoalign) + scalar(@native);
my $file_counter = 0;

if ($total_input == 0) {
	print STDERR sprintf("\n\n");
	print STDERR sprintf("ERROR: No input file specified, you must specify at least 1 bowtie, gsnap, novoalign\n");
	print STDERR sprintf("       or native file for processing.\n");

	print STDERR sprintf("\n");
	exit();
} # End of if statmenet

print STDERR sprintf("\n");		# empty line for formatting purposes

# Remove last slash from temp directory if any
$temp =~ s/\/$//g;

my $run_id = &generateID();
my $total_seconds = 0;
my %chrHandlers;
my ($ifh, $ofh);	# Input and output file handlers
my ($total_acceptable, $total_nr, $total_snps) = (0, 0, 0);

print STDERR sprintf("  o Processing alignment files (%s total)\n", &formatNumber($total_input));

# BOWTIE OUTPUT CONVERSION TO NATIVE SCRIPT FORMAT
for (my $ii=0; $ii < scalar(@bowtie); $ii++) {
	my $starttime = timelocal(localtime(time));		# Current time
	my $lines = 0;
	my $good_reads = 0;
	my $progress = "Please Wait ...";
	
	print STDERR sprintf("     o %s/%s bowtie file : %s - %s", &formatNumber(++$file_counter), &formatNumber($total_input), $bowtie[$ii], $progress);
	
	$ifh = new FileHandle();
	open ($ifh, $bowtie[$ii]) or die("Error opening file\n\n");
	while (<$ifh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my ($name, $strand, $chr, $start, $seq, $qual, $instances, $poly) = split(/\t/, $_);
			
			if ($instances == 0) {
				my $end = $start + length($seq) - 1;

				my @coords = &getBowtiePolymorphisms($start, $end, $strand, $poly);
				my $mis_allowed = scalar(@mismatches) == 1 ? $mismatches[0] : ceil( (length($seq) * $mismatches[0]) / $mismatches[1]);

				if (scalar(@coords) <= $mis_allowed && !exists $unique_reads{$name}) {	# Valid read
					if (!exists $chrHandlers{$chr}) {	# Open a new file handler for this chromosome
						my $tmp = new FileHandle();
						open ($tmp, sprintf(">%s/%s.%s", $temp, $run_id, $chr)) or die("Cannot create temporary file\n");

						# Print headers
						&printTempHeader($tmp);


						$chrHandlers{$chr}->{"handle"} = $tmp;
						$chrHandlers{$chr}->{"filename"} = sprintf("%s/%s.%s", $temp, $run_id, $chr);
						$chrHandlers{$chr}->{"reads"} = 0;
					} # End of if statement

					$ofh = $chrHandlers{$chr}->{"handle"};

					# if assignQual is specified, overwrite existing read qualities
					if (defined($assignQual)) {
						my @qual;
						for (my $i=0; $i < length($seq); $i++) {
							push(@qual, $assignQual);
						} # End of for each statement

						$qual = &phred_to_fastq_sanger(@qual);
					} # End of if statement

					print $ofh sprintf("%s\t%s\t%s\t%s..%s\t%s..%s\t%s\t%s\t%s\t%s\n", $name, "bowtie", 
									   $chr, $start, $end, 1, length($qual), $strand, 
									   $qual, &countNucleotides($seq, "N", "n"), join(" ", @coords));
					$unique_reads{$name} = 1;

					# increment lines count in temp chromsome file
					$chrHandlers{$chr}->{"reads"}++;

					
					$good_reads++;	# Increment good reads counter
				} # End of if statement
			} # End of if statement
		} # End of if statement
	
		if (++$lines % LIMIT == 0) {	# Print lines progress
			$progress = &printProgress($progress, sprintf("%s lines processed", &formatNumber($lines)));
		} # End of if statement
	} # End of while loop
	close ($ifh);
	my $endtime = timelocal(localtime(time));

	my $diff = $endtime - $starttime;
	$total_seconds += $diff;		# Accumulate total run-time
	$progress = &printProgress($progress, sprintf("%s [%s acceptable reads]\n", &formatTime($diff), &formatNumber($good_reads)));

	$total_acceptable += $good_reads;
} # End of for loop

# GSNAP OUTPUT CONVERSION TO NATIVE SCRIPT FORMAT
my ($qualFH, $qualPTR);		# qualFH is the pointer to the quality file, qualPTR is the hash reference for random file access
if (defined($gsnapQuality)) {
	$qualFH = new FileHandle();
	open ($qualFH, $gsnapQuality) or die("cannot open GSNAP quality file $gsnapQuality");
} # end of if statement

for (my $ii=0; $ii < scalar(@gsnap); $ii++) {
	my $starttime = timelocal(localtime(time));		# Current time
	my $reads_count = 0;
	my $good_reads = 0;
	my $progress = "Please Wait ...";

	print STDERR sprintf("     o %s/%s gsnap file : %s - %s", &formatNumber(++$file_counter), &formatNumber($total_input), $gsnap[$ii], $progress);
	
	my ($name, $seq, $qual);
	my ($paths, $paths_num);
	$ifh = new FileHandle();
	open ($ifh, $gsnap[$ii]) or die("Error opening file\n\n");
	while (<$ifh>) {
		chomp;
		if (length($_) != 0) {
			if ($_ =~ m/^>/) {	# Start new entry
				my @fields = split(/\t/, $');

				if (defined($name) && $paths_num != 0) {	# process read
					my $mis_allowed = scalar(@mismatches) == 1 ? $mismatches[0] : ceil( (length($seq) * $mismatches[0]) / $mismatches[1]);
					my $good_path = &evaluateGSNAPPaths($paths, $paths_num, $mis_allowed, $maxSplice);
				
					if (defined($good_path) && !exists $unique_reads{$name}) {
						# compute mismatches and fix EST coordinates to be always from 5' -> 3' according to reference
						my ($genomic_coords, $est_coords, $strand, $poly_coords) = &fixGSNAPPath($name, $seq, $good_path, $paths);
						
						# if assignQual is specified, overwrite existing read qualities and create new
						# quality string using the value in assignQual
						my @quality;
						if (defined($assignQual)) {
							for (my $i=0; $i < length($seq); $i++) {
								push(@quality, $assignQual);
							} # End of for each statement
						
							$qual = &phred_to_fastq_sanger(@quality);
						} # End of if statmenet

						else {	# Fetch quality from file
							if (!defined($qual) && defined($qualFH)) {
								($qualPTR, @quality) = &fetchGSNAPQuality($name, $qualPTR, $qualFH);

								@quality = reverse(@quality) if ($strand eq "-");		# Reverse quality so that its 5' -> 3'
								$qual = &phred_to_fastq_sanger(@quality);
							} # # End of if statement
						} # End of else statement
					
						if (!defined($qual)) {	# No quality found
							print STDERR sprintf("\n\n");
							print STDERR sprintf("ERROR: No quality could be fetched for %s\n", $name);
							print STDERR sprintf("\n");
							exit();
						} # End of if statemnet


						my $chr = $paths->{$good_path}->{"GENOMIC"};
						if (!exists $chrHandlers{$chr}) {	# Open a new file handler for this chromosome
							my $tmp = new FileHandle();
							open ($tmp, sprintf(">%s/%s.%s", $temp, $run_id, $chr)) or die("Cannot create temporary file\n");

							# Print headers

							&printTempHeader($tmp);

							$chrHandlers{$chr}->{"handle"} = $tmp;
							$chrHandlers{$chr}->{"filename"} = sprintf("%s/%s.%s", $temp, $run_id, $chr);
							$chrHandlers{$chr}->{"reads"} = 0;
						} # End of if statement


						$ofh = $chrHandlers{$chr}->{"handle"};

						print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $name, "GSNAP", 
				                   $chr, $genomic_coords, $est_coords, $strand, $qual, 
								   &countNucleotides($seq, "N", "n"), $poly_coords);
						$unique_reads{$name} = 1;
					
						# increment lines count in temp chromsome file
						$chrHandlers{$chr}->{"reads"}++;
						
						$good_reads++;
					} # End of if statement
				} # End of if statement

				if (scalar(@fields) == 4) {
					$name = &trim($fields[3]);
					$qual = &trim($fields[2]);
				} # end of if statement
				else {
					$name = &trim($fields[2]);
					undef($qual);
				} # end of if statemetn

				$seq = $fields[0];
				undef($paths);		# Reset paths
				$paths_num = 0;		# Reset paths count

				if (++$reads_count % LIMIT == 0) {

					$progress = &printProgress($progress, sprintf("%s reads processed", &formatNumber($reads_count)));
				} # End of if statement
			} # End of if statement
			elsif ($_ =~ m/^\s(\S+)/i || $_ =~ m/^,(\S+)/i) {
				my ($alignment, $est_coords, $genomic_coords, $attributes, $scores) = split(/\t/, $_);
				my ($est_start, $est_end) = &getGSNAPESTCoords($est_coords);
				my ($ins, $del, $splice, $sub) = &getGSNAPPolymorphisms($est_start, $est_end, $seq, substr($alignment, 1), $attributes);
				my ($translocation_inversion_scramble) = &GSNAPSpliceTranslocationInversionScramble($attributes);
				my ($strand, $genomic, $gen_start, $gen_end) = &getGSNAPGenomicCoords($genomic_coords);

				if ($alignment =~ m/^\s/) {	# Start new path
					# increment path_num
					$paths_num++;

					# Save values
					$paths->{$paths_num}->{"GENOMIC"} = $genomic;
					$paths->{$paths_num}->{"GENOMIC_COORDS"} = sprintf("%d..%d", $gen_start, $gen_end);
					$paths->{$paths_num}->{"STRAND"} = $strand;
					$paths->{$paths_num}->{"EST_COORDS"} = sprintf("%d..%d", $est_start, $est_end);
					$paths->{$paths_num}->{"EST_TAIL_FRONT"} = $est_start - 1;
					$paths->{$paths_num}->{"EST_TAIL_END"} = length($seq) - $est_end;
					$paths->{$paths_num}->{"INSERTIONS"} = $ins;
					$paths->{$paths_num}->{"DELETIONS"} = $del;
					$paths->{$paths_num}->{"SUBSTITUTIONS"} = $sub;
					$paths->{$paths_num}->{"MAXIMUM_SPLICE_DISTANCE"} = $splice;
					$paths->{$paths_num}->{"TRANSLOCATION_INVERSION_SCRAMBLE"} = $translocation_inversion_scramble;
				} # End of if statement
				else {	# This is a continuation of insertion, deletion, or splice alignment
					# Compute maximum allowable splice distance
					if (exists $paths->{$paths_num}->{"MAXIMUM_SPLICE_DISTANCE"} &&
					    $splice != 0) {
						$paths->{$paths_num}->{"MAXIMUM_SPLICE_DISTANCE"} = &max($paths->{$paths_num}->{"MAXIMUM_SPLICE_DISTANCE"}, $splice);
					} # End of if statemnet

					$paths->{$paths_num}->{"GENOMIC_COORDS"} .= sprintf(" %d..%d", $gen_start, $gen_end);
					$paths->{$paths_num}->{"EST_COORDS"} .= sprintf(" %d..%d", $est_start, $est_end);
					$paths->{$paths_num}->{"SUBSTITUTIONS"} += $sub;	# increment substitutions

					# update tail information
					$paths->{$paths_num}->{"EST_TAIL_END"} = length($seq) - $est_end;
				} # End of else statement

				# Accumulate genomic sequence
				if (exists $paths->{$paths_num}->{"GENOMIC_SEQUENCES"}) {
					$paths->{$paths_num}->{"GENOMIC_SEQUENCES"} .= sprintf(" %s", substr($alignment, 1));
				} # End of if statement
				else {
					$paths->{$paths_num}->{"GENOMIC_SEQUENCES"} = substr($alignment, 1);

				} # End of else statement
			} # End of else if statmeent
		} # End of if statement
	} # End of while loop
	
	# the very last sequence
	if (defined($name) && $paths_num != 0) {	# process read

		my $mis_allowed = scalar(@mismatches) == 1 ? $mismatches[0] : ceil( (length($seq) * $mismatches[0]) / $mismatches[1]);
		my $good_path = &evaluateGSNAPPaths($paths, $paths_num, $mis_allowed, $maxSplice);
	
		if (defined($good_path) && !exists $unique_reads{$name}) {
			# compute mismatches and fix EST coordinates to be always from 5' -> 3' according to reference
			my ($genomic_coords, $est_coords, $strand, $poly_coords) = &fixGSNAPPath($name, $seq, $good_path, $paths);
			
			# if assignQual is specified, overwrite existing read qualities and create new
			# quality string using the value in assignQual

			my @quality;
			if (defined($assignQual)) {
				for (my $i=0; $i < length($seq); $i++) {
					push(@quality, $assignQual);
				} # End of for each statement
			
				$qual = &phred_to_fastq_sanger(@quality);
			} # End of if statmenet
			else {	# Fetch quality from file
				if (!defined($qual) && defined($qualFH)) {
					($qualPTR, @quality) = &fetchGSNAPQuality($name, $qualPTR, $qualFH);
					@quality = reverse(@quality) if ($strand eq "-");		# Reverse quality so that its 5' -> 3'

					$qual = &phred_to_fastq_sanger(@quality);
				} # # End of if statement
			} # End of else statement
		
			if (!defined($qual)) {	# No quality found
				print STDERR sprintf("\n\n");
				print STDERR sprintf("ERROR: No quality could be fetched for %s.\n", $name);
				print STDERR sprintf("\n");
				exit();
			} # End of if statemnet

			my $chr = $paths->{$good_path}->{"GENOMIC"};

			if (!exists $chrHandlers{$chr}) {	# Open a new file handler for this chromosome
				my $tmp = new FileHandle();
				open ($tmp, sprintf(">%s/%s.%s", $temp, $run_id, $chr)) or die("Cannot create temporary file\n");

				# Print headers
				&printTempHeader($tmp);

				$chrHandlers{$chr}->{"handle"} = $tmp;

				$chrHandlers{$chr}->{"filename"} = sprintf("%s/%s.%s", $temp, $run_id, $chr);
				$chrHandlers{$chr}->{"reads"} = 0;
			} # End of if statement

			$ofh = $chrHandlers{$chr}->{"handle"};

			print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $name, "GSNAP", 
					   $chr, $genomic_coords, $est_coords, $strand, $qual, 
					   &countNucleotides($seq, "N", "n"), $poly_coords);
			$unique_reads{$name} = 1;
		
			# increment lines count in temp chromsome file

			$chrHandlers{$chr}->{"reads"}++;
			
			$good_reads++;
		} # End of if statement
	} # End of if statement

	close ($ifh);
	my $endtime = timelocal(localtime(time));
	
	my $diff = $endtime - $starttime;
	$total_seconds += $diff;		# Accumualte total run-time
	$progress = &printProgress($progress, sprintf("%s [%s acceptable reads]\n", &formatTime($diff), &formatNumber($good_reads)));

	$total_acceptable += $good_reads;
} # End of for loop
close ($qualFH) if (defined($qualFH));
undef($qualPTR);	# Free up resources

# NOVOALIGN OUTPUT CONVERSION TO NATIVE SCRIPT FORMAT
for (my $ii=0; $ii < scalar(@novoalign); $ii++) {
	my $starttime = timelocal(localtime(time));		# Current time
	my $lines = 0;
	my $good_reads = 0;
	my $progress = "Please Wait ...";
	
	print STDERR sprintf("     o %s/%s novoalign file : %s - %s", &formatNumber(++$file_counter), &formatNumber($total_input), $novoalign[$ii], $progress);
	
	$ifh = new FileHandle();
	open ($ifh, $novoalign[$ii]) or die("Error opening file\n\n");

	while (<$ifh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my @fields = split(/\t/, $_);

			if ($fields[4] =~ m/^U$/) {	# Only unique hits
				my $name = $fields[0] =~ m/^@/ ? substr($fields[0], 1) : $fields[0];
				my $chr = substr($fields[7], 1);	# reference chromsome
				my $start = $fields[8];
				my $end = $start + length($fields[2]) - 1;
				my $strand = $fields[9] eq "F" ? "+" : "-";
				my $novo_mis = &countNovoalignMismatches($fields[13]);
                my $mis_allowed = scalar(@mismatches) == 1 ? $mismatches[0] : ceil( (length($fields[2]) * $mismatches[0]) / $mismatches[1]);
				my $qual;		# Contains ASCII-33 formatted qualities
				
				if ($novo_mis <= $mis_allowed && !exists $unique_reads{$name}) {
					my @coords = &getNovoalignPolymorphismsCoords($start, $fields[13]);

					if (!exists $chrHandlers{$chr}) {	# Open a new file handler for this chromosome
						my $tmp = new FileHandle();
						open ($tmp, sprintf(">%s/%s.%s", $temp, $run_id, $chr)) or die("Cannot create temporary file\n");

						# Print headers
						&printTempHeader($tmp);


						$chrHandlers{$chr}->{"handle"} = $tmp;
						$chrHandlers{$chr}->{"filename"} = sprintf("%s/%s.%s", $temp, $run_id, $chr);
						$chrHandlers{$chr}->{"reads"} = 0;
					} # End of if statement

					$ofh = $chrHandlers{$chr}->{"handle"};

					# if assignQual is specified, overwrite existing read qualities
					if (defined($assignQual)) {
						my @qual;
						for (my $i=0; $i < length($fields[2]); $i++) {
							push(@qual, $assignQual);
						} # End of for each statement

						$qual = &phred_to_fastq_sanger(@qual);
					} # End of if statement
					else {
						# polymorphism coordinates already from 5' -> 3' from novoalign output,
						# need to flip qualities if orientation is from 3' -> 5'
						$qual = $strand eq "+" ? $fields[3] : join("", reverse(split(//, $fields[3])));
					} # End of else statement

					# Need to adjust reference and read coordinates due to possible
					# insertions/deletions that happens in the read
					my ($ref_coords, $read_coords) = &adjustNovoalignCoordinates($start, $end, 1, length($fields[2]), $fields[13]);

					print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $name, "novoalign", 
									   $chr, $ref_coords, $read_coords, $strand, 
									   $qual, &countNucleotides($fields[2], "N", "n"), join(" ", @coords));

					$unique_reads{$name} = 1;

					# increment lines count in temp chromsome file
					$chrHandlers{$chr}->{"reads"}++;

					
					$good_reads++;	# Increment good reads counter
				} # End of if statement
			} # End of if statement
		} # end of if statement

		if (++$lines % LIMIT == 0) {	# Print lines progress
			$progress = &printProgress($progress, sprintf("%s lines processed", &formatNumber($lines)));
		} # End of if statement
	} # End of while loop
	close ($ifh);
	my $endtime = timelocal(localtime(time));

	my $diff = $endtime - $starttime;
	$total_seconds += $diff;		# Accumulate total run-time
	$progress = &printProgress($progress, sprintf("%s [%s acceptable reads]\n", &formatTime($diff), &formatNumber($good_reads)));

	$total_acceptable += $good_reads;
} # End of for loop

# READING AND SPLITTING NATIVE FILES
for (my $ii=0; $ii < scalar(@native); $ii++) {
	my $starttime = timelocal(localtime(time));		# Current time
	my $lines = 0;
	my $good_reads = 0;
	my $progress = "Please Wait ...";
	
	print STDERR sprintf("     o %s/%s native file : %s - %s", &formatNumber(++$file_counter), &formatNumber($total_input), $native[$ii], $progress);
	
	$ifh = new FileHandle();
	open ($ifh, $native[$ii]) or die("Error opening file\n\n");
	while (<$ifh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my @fields = split(/\t/, $_);
			my $chr = $fields[2];	# reference chromosome

			my $mis_allowed = scalar(@mismatches) == 1 ? $mismatches[0] : ceil( (length($fields[6]) * $mismatches[0]) / $mismatches[1]);
			my $read_mismatches = &getNativeMismatchesCount($fields[$#fields]);

			if (!exists $unique_reads{$fields[0]} && $read_mismatches <= $mis_allowed) {
				if (!exists $chrHandlers{$chr}) {
					my $tmp = new FileHandle();
					open ($tmp, sprintf(">%s/%s.%s", $temp, $run_id, $chr)) or die("Cannot create temporary file\n");

					# Print headers
					&printTempHeader($tmp);
				
					$chrHandlers{$chr}->{"handle"} = $tmp;
					$chrHandlers{$chr}->{"filename"} = sprintf("%s/%s.%s", $temp, $run_id, $chr);
					$chrHandlers{$chr}->{"reads"} = 0;
				} # End of if statement

				$ofh = $chrHandlers{$chr}->{"handle"};

				print $ofh sprintf("%s\n", $_);
				$unique_reads{$fields[0]} = 1;
				
				# increment lines count in temp chromsome file
				$chrHandlers{$chr}->{"reads"}++;
					
				$good_reads++;	# Increment good reads counter
			} # End of if statement
		} # end of if statement

		if (++$lines % LIMIT == 0) {	# Print lines progress
			$progress = &printProgress($progress, sprintf("%s lines processed", &formatNumber($lines)));
		} # End of if statement
	} # End of while loop
	close ($ifh);
	my $endtime = timelocal(localtime(time));

	my $diff = $endtime - $starttime;
	$total_seconds += $diff;		# Accumulate total run-time
	$progress = &printProgress($progress, sprintf("%s [%s acceptable reads]\n", &formatTime($diff), &formatNumber($good_reads)));

	$total_acceptable += $good_reads;
} # End of for loop

print STDERR sprintf("     o TOTAL ACCEPTABLE READS: %s\n", &formatNumber($total_acceptable));

undef(%unique_reads);	# finished reading all output files, don't need read names anymore, free up resources

# OBTAIN ALL CHROMOSOME TO PROCESS
my @chr = sort {$a cmp $b} keys %chrHandlers;

# CLOSE OPEN FILE HANDLERS
foreach my $c (@chr) {
	close ($chrHandlers{$c}->{"handle"});
} # End of for each statement

# READ STACKING FILTERING
if (!$stacking) {
	print STDERR sprintf("  o Reads Stacking Filtering filtering (%s reference sequences)\n", scalar(@chr));

	my ($tlPTR, $logHandle);
	
	if (defined($tlfile)) {
		$logHandle = new FileHandle();
		open ($logHandle, $tlfile) or die("Cannot open trimmed log file\n");
	} # End of if statement

	for (my $i=0; $i < scalar(@chr); $i++) {
		my $starttime = timelocal(localtime(time));
		my $lines = 0;
		my $good_reads = 0;
		my $total = $chrHandlers{$chr[$i]}->{"reads"};		# Total reads in file
		my $progress = "Please Wait ...";
		my (%nr, @keys);

		print STDERR sprintf("     o %d/%d Reference Sequence: %s - %s", $i+1, scalar(@chr), $chr[$i], $progress);

		$ifh = new FileHandle();
		open ($ifh, $chrHandlers{$chr[$i]}->{"filename"}) or die("Cannot open temporary file for $chr[$i]\n");
		while (<$ifh>) {
			chomp;
			if (length($_) != 0 && $_ !~ m/^#/) {
				my ($name, $program, $ref, $ref_coords, $read_coords, $strand, $qual, $n_count, $mis) = split(/\t/, $_);	# Get current line
				my ($ref_start, $ref_end, $mismatches);


				if (defined($tlfile)) {
					($ref_start, $ref_end, $mismatches, $tlPTR) = &getReadAttributesBeforeTrimming($name, $ref_coords, $strand, $mis, $tlPTR, $logHandle);	
				} # end of if statemnet
				else {
					($ref_start, $ref_end, $mismatches) = &getReadAttributes($ref_coords, $mis);
				} # End of else statement

				my $id = sprintf("%d-%d", $ref_start, $ref_end);

				if (!exists $nr{$id}) {	# no other read with same coordinates, save entry
					$nr{$id} = $_;
					push(@keys, $id);	# Save key
					$good_reads++;	# increment by 1

				} # End of if statement
				else {
					# Get saved line
					my ($p_name, $p_program, $p_ref, $p_ref_coords, $p_read_coords, $p_strand, $p_qual, $p_n_count, $p_mis) = split(/\t/, $nr{$id});
					my $p_mismatches;

					if (defined($p_mis)) {
						my @tmp = split(/\s/, $p_mis);
						$p_mismatches = scalar(@tmp);
					} # End of if statement
					else {
						$p_mismatches = 0;
					} # End of else statement

					# validate best non-redundant read, pick longest trimmed sequence with least number of mismatches and least
					# number of n characters
					if (length($qual) > length($p_qual) || (length($qual) == length($p_qual) && $mismatches < $p_mismatches) ||
					    (length($qual) == length($p_qual) && $mismatches == $p_mismatches && $n_count < $p_n_count)) {
						$nr{$id} = $_;	# Replace existing entry with a better one
					} # End of if statement
				} # End of else statemnet

				if (++$lines % LIMIT == 0) {	# Need to print progress
					$progress = &printProgress($progress, sprintf("%s (%2.2f%%) reads filtered", &formatNumber($lines),
											   ($lines / $total) * 100));
				} # End of if statement
			} # End of if statement
		} # End of while loop
		close ($ifh);

		# Writing non-redundant reads
		$progress = &printProgress($progress, "Writing Non-Stackign Reads");

		$ofh = new FileHandle();
		open ($ofh, sprintf(">%s", $chrHandlers{$chr[$i]}->{"filename"})) or die("Cannot write NR reads file\n");
		
		# print headers
		&printTempHeader($ofh);

		foreach my $k (@keys) {
			print $ofh sprintf("%s\n", $nr{$k});
		} # End of for each statemnet
		close ($ofh);
		
		# update non-redundant reads counts
		$chrHandlers{$chr[$i]}->{"NON-STACKING"} = $good_reads;

		my $endtime = timelocal(localtime(time));

		my $diff = $endtime - $starttime;
		$total_seconds += $diff;		# Accumulate total run-time
		$progress = &printProgress($progress, sprintf("%s [%s (%2.2f%%) non-stacking reads remaining]\n", &formatTime($diff),
													  &formatNumber($good_reads), ($good_reads / $total) * 100));
		$total_nr += $good_reads;
	} # End of for loop
	
	print STDERR sprintf("     o TOTAL NON-STACKED READS REMAINING: %s (%2.2f%%)\n", &formatNumber($total_nr), ($total_nr / $total_acceptable) * 100);

	if (defined($logHandle)) {
		close ($logHandle);
		undef($tlPTR);	# free up resources
	} # End of if statement
} # End of if statement

# SNP DISCOVERING STEP
print STDERR sprintf("  o SNP Discovery Step (%s reference sequences)\n", scalar(@chr));

# Print output file headers
$ofh = new FileHandle();
open ($ofh, sprintf(">%s", $output)) or die("Cannot create output file\n");

print $ofh sprintf("##gff-version\t3\n");
print $ofh sprintf("# %s\n", scalar(localtime(time)));
print $ofh sprintf("#\n");
print $ofh sprintf("# SCRIPT NAME: %s\n", $0);
print $ofh sprintf("#\n");
print $ofh sprintf("# INPUT FILES (%d):\n", scalar(@bowtie) + scalar(@gsnap) + scalar(@novoalign) + scalar(@native));

# print bowtie input files
for (my $i=0; $i < scalar(@bowtie); $i++) {
	print $ofh sprintf("#    o %s\n", $bowtie[$i]);
} # End of for loop

# print gsnap input files
for (my $i=0; $i < scalar(@gsnap); $i++) {
	print $ofh sprintf("#    o %s\n", $gsnap[$i]);
} # End of for loop

# print novoalign input files
for (my $i=0; $i < scalar(@novoalign); $i++) {
	print $ofh sprintf("#    o %s\n", $novoalign[$i]);
} # End of for loop

# print native input files
for (my $i=0; $i < scalar(@native); $i++) {
	print $ofh sprintf("#    o %s\n", $native[$i]);
} # End of for loop

print $ofh sprintf("#\n");

if (defined($gsnapQuality)) {
	print $ofh sprintf("# GSNAP QUALITY FILE: %s\n", $gsnapQuality);

} # End of if statmeent

if (defined($assignQual)) {
	print $ofh sprintf("# DEFAULT QUALITY VALUE: %s\n", $assignQual);
} # End of if statmenet

if (scalar(@mismatches) == 2) {
	print $ofh sprintf("# MAXIMUM MISMATCHES ALLOWED PER READ: %s %s (allow at most %s mismatches per %s bp)\n",
	                   $mismatches[0], $mismatches[1], $mismatches[0], $mismatches[1]);
} # End of if statement
else {
	print $ofh sprintf("# MAXIMUM MISMATCHES ALLOWED PER READ: %s (regardless of read length)\n",
	                   $mismatches[0]);
} # End of else statement

print $ofh sprintf("# MAXIMUM GSNAP SPLICE DISTANCE ALLOWED: %s\n", $maxSplice) if (scalar(@gsnap) != 0);
print $ofh sprintf("# MAXIMUM NUMBER OF ALLELES PER SITE: %s\n", $total_alleles);
print $ofh sprintf("# MINIMUM READ COUNT ALLELE: %d\n", $num);
print $ofh sprintf("# MINIMUM READ COVERAGE PER ALLELE: %2.2f\n", $allele_coverage);
print $ofh sprintf("# MINIMUM OVERALL ALLELE COVERAGE: %2.2f\n", $overall_coverage);
print $ofh sprintf("# MINIMUM ALLELE QUALITY VALUE: %s\n", $min_quality);
print $ofh sprintf("# BASES IGNORED AT BEGINNING/END OF READ: %s\n", $ignore);
print $ofh sprintf("# FILTER STACKING READS: %s\n", $stacking ? "NO" : "YES");
print $ofh sprintf("# TRIMMED LOG FILE: %s\n", $tlfile) if (defined($tlfile) && !$stacking);
print $ofh sprintf("#\n");
print $ofh sprintf("# NOTE: SNP Discovery statistics located at the end of the file\n");
print $ofh sprintf("\n");

for (my $ii=0; $ii < scalar(@chr); $ii++) {
	my $starttime = timelocal(localtime(time));
	my $lines = 0;
	my $total = !$stacking ? $chrHandlers{$chr[$ii]}->{"NON-STACKING"} : $chrHandlers{$chr[$ii]}->{"reads"};
	my $progress = "Please Wait ...";
	my (%snp, %pm, %ref, @sites);	# variables for SNPs, perfect match counts, reference alleles, snp sites
	
	print STDERR sprintf("     o %d/%d Talling Polymorphisms: %s - %s", $ii+1, scalar(@chr), $chr[$ii], $progress);

	$ifh = new FileHandle();
	open ($ifh, $chrHandlers{$chr[$ii]}->{"filename"}) or die("cannot open temporary chromosome file\n");
	
	while (<$ifh>) {
		chomp;
		if (length($_) != 0 && $_ !~ m/^#/) {
			my ($name, $program, $ref, $ref_coords, $read_coords, $strand, $qual, $n_count, $mis) = split(/\t/, $_);	# Get current line
			my @ref_pieces = split(/\s/, $ref_coords);
			my @read_pieces = split(/\s/, $read_coords);

			my @mismatches = split(/\s/, $mis);
			my @qualities = &fastq_sanger_to_phred($qual);	# convert qualities to ASCII-33 format
			my ($start, $end);

			# Talling perfect match nucleotides
			for (my $i=0; $i < scalar(@ref_pieces); $i++) {
				my ($ref_piece_start, $ref_piece_end) = split(/\.\./, $ref_pieces[$i]);
				my ($read_piece_start, $read_piece_end) = split(/\.\./, $read_pieces[$i]);

				# adjusting bases to be ignored at the beginning of read
				if ($i == 0) {
					$ref_piece_start += $ignore;

					$read_piece_start += $ignore;

				} # End of if statement

				# adjusting bases to be ignored at the end of the read
				if ($i == $#ref_pieces) {
					$ref_piece_end -= $ignore;
					$read_piece_end -= $ignore;
				} # End of if statement

				# Consider SNPs in this start/end range only
				$start = defined($start) ? &min($ref_piece_start, $start) : $ref_piece_start;
				$end = defined($end) ? &max($ref_piece_end, $end) : $ref_piece_end;

				my $offset = 0;
				for (my $j=$ref_piece_start; $j <= $ref_piece_end; $j++) {
					my $adj_read_pos = $read_piece_start + $offset;
					if ($qualities[$adj_read_pos - 1] >= $min_quality) {
						if (exists $pm{$j}) {
							$pm{$j}++;
						} # end of if statement

						else {
							$pm{$j} = 1;

						} # End of else statement
					} # End of if statement
					$offset++;
				} # End of for loop
			} # End of for each statement
	
			# Talling mismatches
			foreach my $m (@mismatches) {
				if ($m =~ m/^(\d+):(\S)>(\d+):(\S)$/) {
					my $ref_pos = $1;
					my $ref_allele = $2;
					my $read_pos = $3;
					my $read_allele = $4;

					if ($qualities[$read_pos - 1] >= $min_quality && $ref_pos >= $start && $ref_pos <= $end) {
						# decrement perfect match count
						$pm{$ref_pos}--;

						# tally polymorphism
						if (!exists $snp{$ref_pos}) {
							$ref{$ref_pos} = $ref_allele;
							$snp{$ref_pos} = $read_allele;
							push(@sites, $ref_pos);
						} # End of if statement

						else {
							$snp{$ref_pos} .= $read_allele;
						} # end of else statemen

					} # End of if statement
				} # End of if statement
				else {
					print STDERR sprintf("\n\n");
					print STDERR sprintf("ERROR: THIS SHOULDN'T HAPPEN\n");
					print STDERR sprintf("\n");
					exit();
				} # End of else statement
			} # End of for each statement
		
			if (++$lines % LIMIT == 0) {	# Need to print progress
				$progress = &printProgress($progress, sprintf("%s (%2.2f%%) reads processed", &formatNumber($lines),
										   ($lines / $total) * 100));
			} # End of if statement
		} # End of if statement
	} # End of while loop
	close ($ifh);

	# Processing each possible SNP sites
	$progress = &printProgress($progress, sprintf("Processing %s potential SNP Sites. Please Wait ...", &formatNumber(scalar(@sites))));

	my $good_snps = 0;
	foreach my $pos (sort {$a <=> $b} @sites) {
		if (exists $ref{$pos}) {
			my $ref_allele = uc($ref{$pos});
			my $read_alleles = $snp{$pos};
			my $perfect_match = exists $pm{$pos} ? $pm{$pos} : 0;

			for (my $i=0; $i < $perfect_match; $i++) {
				$read_alleles .= $ref_allele;
			} # End of for loop

			$read_alleles = uc($read_alleles);
			my @rank = &rankAlleles($read_alleles, $total_alleles);
			my @allele_groups;
			my $total_allelic_reads = 0;

			foreach my $ra (@rank) {
				if (length($ra) >= $num && (length($ra) / length($read_alleles)) >= $allele_coverage) {
					push(@allele_groups, $ra);
					$total_allelic_reads += length($ra);
				} # End of if statemnet
			} # End of for each statement

			if ($total_allelic_reads / length($read_alleles) >= $overall_coverage && &countNonReferenceAlelles($ref_allele, @allele_groups) >= 1) {
				print $ofh sprintf("%s\t%s\t%s\t%s\t%s\t.\t.\t.\tName=%s-%s_%d;Note=%s",
				                   $chr[$ii], $source, $feature, $pos, $pos, $source, $chr[$ii], $pos, $ref_allele);

				# printing SNP alleles
				foreach my $ag (@allele_groups) {
					my $nt = substr($ag, 0, 1);
					print $ofh sprintf("/%s", substr($ag, 0, 1)) if ($nt ne $ref_allele);
				} # End of for each statement

				# printing allele gorups info
				print $ofh sprintf(", Ref: %s, Alleles: %d [", $ref_allele, scalar(@allele_groups));

				for (my $i=0; $i < scalar(@allele_groups); $i++) {
					my $nt = substr($allele_groups[$i], 0, 1);
					my $nt_count = length($allele_groups[$i]);
					my $cvg = ($nt_count / length($read_alleles)) * 100;

					if ($i == 0) {
						print $ofh sprintf("%s (%s reads - %2.2f%%)", $nt, $nt_count, $cvg);
					} # End of if statement
					else {
						print $ofh sprintf(", %s (%s reads - %2.2f%%)", $nt, $nt_count, $cvg);
					} # end of else statement
				} # End of for each statement

				print $ofh sprintf("], %s reads total\n", length($read_alleles));
				$good_snps++;
			} # End of if statement
		} # End of if statement
	} # End of for each statement

	# update SNP counts
	$chrHandlers{$chr[$ii]}->{"SNPs"} = $good_snps;

	my $endtime = timelocal(localtime(time));

	my $diff = $endtime - $starttime;
	$total_seconds += $diff;		# Accumulate total run-time
	$progress = &printProgress($progress, sprintf("%s [%s SNPs Discovered]\n", &formatTime($diff), &formatNumber($good_snps)));
	$total_snps += $good_snps;
} # End of for loop

# Writing SNP discovery statistics to output file
print $ofh sprintf("\n");
print $ofh sprintf("# TOTAL ACCEPTABLE READS: %s\n", &formatNumber($total_acceptable));
print $ofh sprintf("# TOTAL NON-REDUNDANT READS: %s (%2.2f%%)\n", &formatNumber($total_nr), ($total_nr / $total_acceptable) * 100) if (!$stacking);
print $ofh sprintf("# SNPs DISCOVERED PER REFERENCE:\n");
for (my $i=0; $i < scalar(@chr); $i++) {
	print $ofh sprintf("#     o %s - %s (%2.2f%%)\n", $chr[$i], &formatNumber($chrHandlers{$chr[$i]}->{"SNPs"}),
	                    ($chrHandlers{$chr[$i]}->{"SNPs"} / $total_snps) * 100);
} # End of for statement
print $ofh sprintf("# SNPs DISCOVERED: %s\n", &formatNumber($total_snps));
print $ofh sprintf("# TOTAL RUN-TIME: %s\n", &formatTime($total_seconds));
close ($ofh);

print STDERR sprintf("  o TOTAL SNPs DISCOVERED: %s\n", &formatNumber($total_snps));
print STDERR sprintf("  o TOTAL RUN-TIME: %s\n", &formatTime($total_seconds));
print STDERR sprintf("\n");

# Clean up temporary files
if ($clean) {
	foreach my $c (@chr) {
		unlink($chrHandlers{$c}->{"filename"});
	} # End of for each statement
} # End of for each statement
