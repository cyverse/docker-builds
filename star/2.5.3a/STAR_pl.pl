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

use English        qw(-no_match_vars);
use File::Basename qw(basename);
use File::Copy     qw(move);
use File::Path     qw(make_path);
use File::Spec;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Readonly;

our $VERSION = '1.0.0';

Readonly my $THREADS   => 4;
Readonly my $INDEX_DIR => 'index';
Readonly my $OUT_DIR   => 'output';
Readonly my $BAM_DIR   => 'bam_output';

# Field layout of a wait status, as returned in $CHILD_ERROR; see perlvar.
Readonly my $EXEC_FAILED     => -1;
Readonly my $SIGNAL_MASK     => 127;
Readonly my $EXIT_CODE_SHIFT => 8;

exit main();

# Parses the command line and orchestrates the calls to STAR.
sub main {
    my $option = parse_command_line();
    validate_arguments($option);

    my @annotation_args = annotation_args( $option->{user_annotation} );

    build_index(
        {   database        => $option->{user_database},
            annotation_args => \@annotation_args,
        }
    );

    make_path( $OUT_DIR, $BAM_DIR );

    my @queries = @{ $option->{file_query} };
    my @mates   = @{ $option->{file_query2} };
    my $paired  = is_paired_end( $option->{file_type} );

    for my $index ( 0 .. $#queries ) {
        align(
            {   read_file       => $queries[$index],
                mate_file       => $paired ? $mates[$index] : undef,
                annotation_args => \@annotation_args,
                star_args       => $option->{star_args},
            }
        );
    }

    return 0;
}

# Parses the command line.
sub parse_command_line {
    my %option = (
        file_query      => [],
        file_query2     => [],
        user_database   => undef,
        user_annotation => undef,
        file_type       => undef,
    );

    GetOptions(
        'file_query=s'      => $option{file_query},
        'file_query2=s'     => $option{file_query2},
        'user_database=s'   => \$option{user_database},
        'user_annotation=s' => \$option{user_annotation},
        'file_type=s'       => \$option{file_type},
    ) or die "Error: unable to parse the command line\n";

    # Whatever Getopt::Long left behind is handed straight to STAR.
    $option{star_args} = [@ARGV];

    return \%option;
}

# Validates the command-line options in order to provide useful error messages
# in cases where the analysis would fail.
sub validate_arguments {
    my ($option) = @_;

    my @queries = @{ $option->{file_query} };
    my @mates   = @{ $option->{file_query2} };

    if ( !@queries ) {
        die "Error: no FASTQ files were supplied\n";
    }

    my $database = $option->{user_database};
    if ( !defined $database || $database eq q{} ) {
        die "Error: no reference genome was supplied\n";
    }
    if ( !-f $database ) {
        die "Error: the reference genome $database does not exist\n";
    }
    if ( !looks_like_fasta($database) ) {
        die "Error: the reference genome $database is not a FASTA file\n";
    }

    my $annotation = $option->{user_annotation};
    if ( defined $annotation && !-f $annotation ) {
        die "Error: the annotation $annotation does not exist\n";
    }

    for my $query_file ( @queries, @mates ) {
        if ( !-f $query_file ) {
            die "Error: the FASTQ file $query_file does not exist\n";
        }
    }

    validate_file_type( $option->{file_type}, \@queries, \@mates );

    return;
}

# Verifies that the file type argument is recognized and compatible with the
# input files that were provided.
sub validate_file_type {
    my ( $file_type, $queries, $mates ) = @_;

    if ( !defined $file_type ) {
        die "Error: no file type was supplied; expected SE or PE\n";
    }

    my $type = uc $file_type;
    if ( $type ne 'SE' && $type ne 'PE' ) {
        die "Error: unrecognized file type '$file_type'; expected SE or PE\n";
    }

    if ( $type eq 'SE' ) {
        if ( @{$mates} ) {
            die 'Error: --file_query2 was supplied but the file type is SE; '
                . "use --file_type PE to align these files as paired ends\n";
        }

        return;
    }

    if ( !@{$mates} ) {
        die "Error: at least one file for each paired end is required\n";
    }
    if ( @{$queries} != @{$mates} ) {
        die "Error: unequal number of files for paired ends\n";
    }

    return;
}

# Returns true if the user requested a paired-end alignment.
sub is_paired_end {
    my ($file_type) = @_;

    return ( uc $file_type ) eq 'PE';
}

# Returns true if the file appears to be a FASTA file.
sub looks_like_fasta {
    my ($path) = @_;

    my $first_line = q{};

    open my $fh, '<', $path
        or die "Error: cannot read $path: $OS_ERROR\n";
    while ( my $line = <$fh> ) {
        next if $line =~ /\A \s* \z/xms;
        $first_line = $line;
        last;
    }
    close $fh
        or die "Error: cannot close $path: $OS_ERROR\n";

    return $first_line =~ /\A > /xms ? 1 : 0;
}

# Returns the STAR command-line option for the annotation if the user provided
# an annotation file.
sub annotation_args {
    my ($annotation) = @_;

    if ( !defined $annotation ) {
        return;
    }

    return ( '--sjdbGTFfile', $annotation );
}

# Calls STAR in order to build the index.
sub build_index {
    my ($arg) = @_;

    my $database        = $arg->{database};
    my $annotation_args = $arg->{annotation_args};

    my $name = basename( $database, qw(.fa .fas .fasta .fna) );
    report("STAR-indexing $name");

    make_path($INDEX_DIR);

    run('STAR',
        '--runThreadN'       => $THREADS,
        '--runMode'          => 'genomeGenerate',
        '--genomeDir'        => $INDEX_DIR,
        '--genomeFastaFiles' => $database,
        @{$annotation_args},
    );

    my $state
        = @{$annotation_args} ? 'with_annotation' : 'without_annotation';
    print {*STDOUT} "index $state\n"
        or die "Error: cannot write to standard output: $OS_ERROR\n";

    return;
}

# Calls STAR in order to perform the alignment then moves output files into
# per-sample and alignment destination directories.
sub align {
    my ($arg) = @_;

    my $read_file = $arg->{read_file};
    my $mate_file = $arg->{mate_file};

    my $prefix = basename($read_file);
    $prefix =~ s/[.] \S+ \z//xms;

    my @read_files = ($read_file);
    if ( defined $mate_file ) {
        push @read_files, $mate_file;
    }

    my @align_command = (
        'STAR',
        @{ $arg->{annotation_args} },
        @{ $arg->{star_args} },
        '--runThreadN'        => $THREADS,
        '--genomeDir'         => $INDEX_DIR,
        '--outReadsUnmapped'  => 'Fastx',
        '--outFileNamePrefix' => "$prefix.",
        '--readFilesIn'       => @read_files,
        '--readFilesCommand'  => 'gunzip',
        '-c',
    );

    report("Executing: @align_command");
    run(@align_command);

    collect_results($prefix);

    return;
}

# Moves output files from STAR into per-sample and alignment directories.
sub collect_results {
    my ($prefix) = @_;

    my $sample_dir = File::Spec->catdir( $OUT_DIR, $prefix );
    make_path($sample_dir);

    opendir my $dh, q{.}
        or die "Error: cannot read the working directory: $OS_ERROR\n";
    my @entries = grep { $_ ne q{.} && $_ ne q{..} } readdir $dh;
    closedir $dh
        or die "Error: cannot close the working directory: $OS_ERROR\n";

    my %is_reserved = map { $_ => 1 } ( $INDEX_DIR, $OUT_DIR, $BAM_DIR );

    for my $entry (@entries) {
        next if $is_reserved{$entry};

        my $destination = destination_for( $entry, $prefix, $sample_dir );
        next if !defined $destination;

        move( $entry, File::Spec->catfile( $destination, $entry ) )
            or warn
            "Warning: could not move $entry into $destination: $OS_ERROR\n";
    }

    return;
}

# Returns the destination for a file or returns `undef` if the file shouldn't
# be moved.
sub destination_for {
    my ( $entry, $prefix, $sample_dir ) = @_;

    if (   $entry =~ /\A Log/xms
        || $entry =~ /STARgenome \z/xms
        || $entry =~ /\A \Q$prefix\E .* out \z/xms
        || $entry =~ /tab \z/xms
        || $entry =~ /Unmapped/xms )
    {
        return $sample_dir;
    }

    return $BAM_DIR if $entry =~ /[.] (?: bam | sam ) \z/xms;

    return;
}

# Runs a command as a subprocess, exiting if the command fails.
sub run {
    my @command = @_;

    if ( system(@command) != 0 ) {
        die "Error: @command failed: " . describe_exit($CHILD_ERROR) . "\n";
    }

    return;
}

# Returns a description of a wait status.
sub describe_exit {
    my ($status) = @_;

    if ( $status == $EXEC_FAILED ) {
        return "could not be executed: $OS_ERROR";
    }

    my $signal = $status & $SIGNAL_MASK;
    if ($signal) {
        return "killed by signal $signal";
    }

    return 'exit status ' . ( $status >> $EXIT_CODE_SHIFT );
}

# Prints a message to stderr, exiting if the write fails.
sub report {
    my ($message) = @_;

    print {*STDERR} "$message\n"
        or die "Error: cannot write to standard error: $OS_ERROR\n";

    return;
}
