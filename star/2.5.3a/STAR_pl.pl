#!/usr/bin/env perl
#
# STAR_pl.pl -- wrapper around STAR 2.5.3a.
#
# Builds a genome index from a user-supplied FASTA (optionally with a GTF
# annotation), then aligns each supplied FASTQ, single- or paired-end.
# Per-sample logs and auxiliary files are collected under output/<basename>/
# and the alignments under bam_output/.
#
# Any arguments not consumed below are passed through to STAR verbatim.

use strict;
use warnings;

use File::Basename qw(basename);
use File::Copy qw(move);
use File::Path qw(make_path);
use File::Spec;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use constant {
    THREADS   => 4,
    INDEX_DIR => 'index',
    OUT_DIR   => 'output',
    BAM_DIR   => 'bam_output',
};

my (@file_query, @file_query2, $user_database_path, $user_annotation_path, $file_type);

GetOptions(
    "file_query=s"      => \@file_query,
    "file_query2=s"     => \@file_query2,
    "user_database=s"   => \$user_database_path,
    "user_annotation=s" => \$user_annotation_path,
    "file_type=s"       => \$file_type,
) or die "Error: unable to parse the command line\n";

# Whatever Getopt::Long left behind is handed straight to STAR.
my @star_args = @ARGV;

my $format = validate_arguments();

my @annotation_args = defined $user_annotation_path
    ? ('--sjdbGTFfile', $user_annotation_path)
    : ();

build_index();

make_path(OUT_DIR, BAM_DIR);

for my $i (0 .. $#file_query) {
    align($file_query[$i], ($format eq 'PE' ? $file_query2[$i] : undef));
}

exit 0;

# Checks every input up front and returns the normalized file type, so that a
# bad invocation fails before STAR spends an hour building an index.
sub validate_arguments {
    @file_query
        or die "Error: no FASTQ files were supplied\n";

    defined $user_database_path && length $user_database_path
        or die "Error: no reference genome was supplied\n";

    -f $user_database_path
        or die "Error: the reference genome $user_database_path does not exist\n";

    looks_like_fasta($user_database_path)
        or die "Error: the reference genome $user_database_path is not a FASTA file\n";

    if (defined $user_annotation_path) {
        -f $user_annotation_path
            or die "Error: the annotation $user_annotation_path does not exist\n";
    }

    for my $query_file (@file_query, @file_query2) {
        -f $query_file
            or die "Error: the FASTQ file $query_file does not exist\n";
    }

    defined $file_type
        or die "Error: no file type was supplied; expected SE or PE\n";

    my $normalized = uc $file_type;
    $normalized eq 'SE' || $normalized eq 'PE'
        or die "Error: unrecognized file type '$file_type'; expected SE or PE\n";

    if ($normalized eq 'PE' || @file_query2) {
        @file_query && @file_query2
            or die "Error: at least one file for each paired end is required\n";
        @file_query == @file_query2
            or die "Error: unequal number of files for paired ends\n";
    }

    return $normalized;
}

# Reads only as far as the first non-blank line rather than shelling out to
# grep, which would interpret the path as part of a shell command.
sub looks_like_fasta {
    my ($path) = @_;

    open my $fh, '<', $path
        or die "Error: cannot read $path: $!\n";

    while (my $line = <$fh>) {
        next if $line =~ /^\s*$/;
        close $fh;
        return $line =~ /^>/;
    }

    close $fh;
    return 0;
}

sub build_index {
    my $name = basename($user_database_path, qw(.fa .fas .fasta .fna));
    report("STAR-indexing $name");

    make_path(INDEX_DIR);

    run(
        'STAR',
        '--runThreadN',        THREADS,
        '--runMode',           'genomeGenerate',
        '--genomeDir',         INDEX_DIR,
        '--genomeFastaFiles',  $user_database_path,
        @annotation_args,
    );

    print defined $user_annotation_path
        ? "index with_annotation\n"
        : "index without_annotation\n";
}

sub align {
    my ($query_file, $second_file) = @_;

    my $basename = basename($query_file);
    $basename =~ s/\.\S+$//;

    my @read_files = ($query_file);
    push @read_files, $second_file if defined $second_file;

    my @align_command = (
        'STAR',
        @annotation_args,
        @star_args,
        '--runThreadN',        THREADS,
        '--genomeDir',         INDEX_DIR,
        '--outReadsUnmapped',  'Fastx',
        '--outFileNamePrefix', "$basename.",
        '--readFilesIn',       @read_files,
        '--readFilesCommand',  'gunzip', '-c',
    );

    report("Executing: @align_command");
    run(@align_command);

    collect_results($basename);
}

# STAR writes its output into the working directory, so sort that directory by
# hand instead of handing shell globs to mv.
sub collect_results {
    my ($basename) = @_;

    my $sample_dir = File::Spec->catdir(OUT_DIR, $basename);
    make_path($sample_dir);

    opendir my $dh, '.'
        or die "Error: cannot read the working directory: $!\n";
    my @entries = grep { $_ ne '.' && $_ ne '..' } readdir $dh;
    closedir $dh;

    my %reserved = map { $_ => 1 } (INDEX_DIR, OUT_DIR, BAM_DIR);

    for my $entry (@entries) {
        next if $reserved{$entry};

        my $destination = destination_for($entry, $basename, $sample_dir);
        next unless defined $destination;

        move($entry, File::Spec->catfile($destination, $entry))
            or warn "Warning: could not move $entry into $destination: $!\n";
    }
}

sub destination_for {
    my ($entry, $basename, $sample_dir) = @_;

    return $sample_dir
        if $entry =~ /^Log/
        || $entry =~ /STARgenome$/
        || $entry =~ /^\Q$basename\E.*out$/
        || $entry =~ /tab$/
        || $entry =~ /Unmapped/;

    return BAM_DIR if $entry =~ /bam$/;

    return;
}

# Runs a command as an argument list, so no part of it is interpreted by a
# shell, and stops on failure instead of pressing on with missing output.
sub run {
    my @command = @_;

    return if system(@command) == 0;

    die "Error: @command failed: " . describe_exit($?) . "\n";
}

sub describe_exit {
    my ($status) = @_;

    return "could not be executed: $!" if $status == -1;
    return 'killed by signal ' . ($status & 127) if $status & 127;
    return 'exit status ' . ($status >> 8);
}

sub report {
    my ($message) = @_;
    print STDERR "$message\n";
}
