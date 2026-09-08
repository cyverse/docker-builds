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

# Owns every value the run depends on and hands each subroutine exactly what
# it needs, so that nothing below reaches outside its own scope for input.
sub main {
    my $option = parse_command_line();
    my $format = validate_arguments($option);

    my @annotation_args = annotation_args( $option->{user_annotation} );

    build_index(
        {   database        => $option->{user_database},
            annotation_args => \@annotation_args,
        }
    );

    make_path( $OUT_DIR, $BAM_DIR );

    my @queries = @{ $option->{file_query} };
    my @mates   = @{ $option->{file_query2} };

    for my $index ( 0 .. $#queries ) {
        align(
            {   read_file       => $queries[$index],
                mate_file       => $format eq 'PE' ? $mates[$index] : undef,
                annotation_args => \@annotation_args,
                star_args       => $option->{star_args},
            }
        );
    }

    return 0;
}

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

# Checks every input up front and returns the normalized file type, so that a
# bad invocation fails before STAR spends an hour building an index.
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

    return normalize_file_type( $option->{file_type}, \@queries, \@mates );
}

sub normalize_file_type {
    my ( $file_type, $queries, $mates ) = @_;

    if ( !defined $file_type ) {
        die "Error: no file type was supplied; expected SE or PE\n";
    }

    my $normalized = uc $file_type;
    if ( $normalized ne 'SE' && $normalized ne 'PE' ) {
        die "Error: unrecognized file type '$file_type'; expected SE or PE\n";
    }

    if ( $normalized eq 'PE' || @{$mates} ) {
        if ( !@{$queries} || !@{$mates} ) {
            die "Error: at least one file for each paired end is required\n";
        }
        if ( @{$queries} != @{$mates} ) {
            die "Error: unequal number of files for paired ends\n";
        }
    }

    return $normalized;
}

# Reads only as far as the first non-blank line rather than shelling out to
# grep, which would interpret the path as part of a shell command.
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

# Returns the STAR flags that name the annotation, or an empty list when the
# run has none, so that callers can interpolate the result unconditionally.
sub annotation_args {
    my ($annotation) = @_;

    if ( !defined $annotation ) {
        return;
    }

    return ( '--sjdbGTFfile', $annotation );
}

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

# STAR writes its output into the working directory, so sort that directory by
# hand instead of handing shell globs to mv.
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

    return $BAM_DIR if $entry =~ /bam \z/xms;

    return;
}

# Runs a command as an argument list, so no part of it is interpreted by a
# shell, and stops on failure instead of pressing on with missing output.
sub run {
    my @command = @_;

    if ( system(@command) != 0 ) {
        die "Error: @command failed: " . describe_exit($CHILD_ERROR) . "\n";
    }

    return;
}

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

sub report {
    my ($message) = @_;

    print {*STDERR} "$message\n"
        or die "Error: cannot write to standard error: $OS_ERROR\n";

    return;
}
