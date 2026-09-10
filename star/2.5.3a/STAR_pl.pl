#!/usr/bin/env perl
#
# STAR_pl.pl -- wrapper around STAR 2.5.3a.
#
# Builds a genome index from a user-supplied FASTA (optionally with a GTF
# annotation), then aligns each supplied FASTQ, single- or paired-end.
# Per-sample logs and auxiliary files are collected under output/<basename>/
# and the alignments under bam_output/.
#
# STAR's memory use is capped at the memory the job has, less a margin for
# everything else running in the container. That figure comes from
# --memory_limit (in GiB) if it is given, and otherwise from the container's
# own memory limit.
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
use List::Util   qw(any min sum0);
use Readonly;

our $VERSION = '1.0.0';

Readonly my $THREADS   => 4;
Readonly my $INDEX_DIR => 'index';
Readonly my $OUT_DIR   => 'output';
Readonly my $BAM_DIR   => 'bam_output';

# Units of memory, and how much of the job's memory is kept back for the
# processes in the container other than STAR.
Readonly my $BYTES_PER_KIB   => 1024;
Readonly my $BYTES_PER_GIB   => $BYTES_PER_KIB**3;
Readonly my $MEMORY_HEADROOM => $BYTES_PER_GIB;

# Where Linux reports the cgroups of a process and the memory of the machine.
Readonly my $PROC_CGROUP  => '/proc/self/cgroup';
Readonly my $PROC_MEMINFO => '/proc/meminfo';
Readonly my $CGROUP_ROOT  => '/sys/fs/cgroup';

# Field layout of a wait status, as returned in $CHILD_ERROR; see perlvar.
Readonly my $EXEC_FAILED     => -1;
Readonly my $SIGNAL_MASK     => 127;
Readonly my $EXIT_CODE_SHIFT => 8;

exit main();

# Parses the command line and orchestrates the calls to STAR.
sub main {
    my $option = parse_command_line();
    validate_arguments($option);

    my $star_memory     = memory_for_star( $option->{memory_limit} );
    my @annotation_args = annotation_args( $option->{user_annotation} );

    build_index(
        {   database        => $option->{user_database},
            annotation_args => \@annotation_args,
            memory          => $star_memory,
        }
    );
    check_index_fits($star_memory);

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
                memory          => $star_memory,
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
        memory_limit    => undef,
    );

    GetOptions(
        'file_query=s'      => $option{file_query},
        'file_query2=s'     => $option{file_query2},
        'user_database=s'   => \$option{user_database},
        'user_annotation=s' => \$option{user_annotation},
        'file_type=s'       => \$option{file_type},
        'memory_limit=s'    => \$option{memory_limit},
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

    my $memory_limit = $option->{memory_limit};
    if ( defined $memory_limit && !is_positive_number($memory_limit) ) {
        die "Error: invalid memory limit '$memory_limit'; "
            . "expected a number of GiB such as 16 or 7.5\n";
    }

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

# Returns true if the value is a plain decimal number greater than zero.
sub is_positive_number {
    my ($value) = @_;

    return $value =~ /\A \d+ (?: [.] \d+ )? \z/xms && $value > 0;
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

# Returns the number of bytes of memory that STAR may use: the memory the job
# has, less the headroom kept for everything else in the container.
sub memory_for_star {
    my ($requested_gib) = @_;

    my ( $memory, $source ) = memory_available($requested_gib);

    my $star_memory = $memory - $MEMORY_HEADROOM;
    if ( $star_memory <= 0 ) {
        die 'Error: the job has '
            . gib($memory)
            . ' of memory, which leaves nothing for STAR after keeping '
            . gib($MEMORY_HEADROOM)
            . " for other processes\n";
    }

    report(   'Memory available to the job: '
            . gib($memory)
            . " ($source); STAR may use "
            . gib($star_memory) );

    return $star_memory;
}

# Returns the memory the job has, in bytes, along with a description of where
# the figure came from.
sub memory_available {
    my ($requested_gib) = @_;

    my $limit = container_memory_limit();

    if ( defined $requested_gib ) {
        my $requested = int( $requested_gib * $BYTES_PER_GIB );
        if ( defined $limit && $requested > $limit ) {
            warn
                "Warning: --memory_limit $requested_gib GiB is more than the "
                . 'container memory limit of '
                . gib($limit)
                . "; using the container limit instead\n";
            return ( $limit, 'container memory limit' );
        }

        return ( $requested, '--memory_limit' );
    }

    if ( defined $limit ) {
        return ( $limit, 'container memory limit' );
    }

    my $available = meminfo_bytes('MemAvailable');
    if ( defined $available ) {
        return ( $available, 'memory currently available on the host' );
    }

    die 'Error: cannot determine how much memory this job has; '
        . "please specify --memory_limit\n";
}

# Returns the memory limit, in bytes, that the cgroups of this process impose,
# or undef if they impose none.
sub container_memory_limit {
    my @limits = map { read_limit($_) } cgroup_limit_files();

    # A limit of at least the size of the machine constrains nothing. This is
    # also how cgroup v1 reports an unlimited cgroup.
    my $total = meminfo_bytes('MemTotal');
    if ( defined $total ) {
        @limits = grep { $_ < $total } @limits;
    }

    if ( !@limits ) {
        return;
    }

    return min(@limits);
}

# Returns the paths of the memory limit files for the cgroups of this process
# and all of their ancestors, under both cgroup v2 and cgroup v1. A parent
# cgroup can impose a tighter limit than the process's own, so every level
# counts, and levels that do not exist in this container are skipped later.
sub cgroup_limit_files {
    open my $fh, '<', $PROC_CGROUP
        or return;
    my @memberships = <$fh>;
    close $fh
        or die "Error: cannot close $PROC_CGROUP: $OS_ERROR\n";

    my @files;
    for my $membership (@memberships) {
        my ( $hierarchy, $controllers, $path )
            = $membership =~ /\A ([^:]*) : ([^:]*) : (.*?) \s* \z/xms;
        next if !defined $path;

        my ( $base, $file );
        if ( $hierarchy eq '0' && $controllers eq q{} ) {
            ( $base, $file ) = ( $CGROUP_ROOT, 'memory.max' );
        }
        elsif ( any { $_ eq 'memory' } split /,/xms, $controllers ) {
            ( $base, $file )
                = ( "$CGROUP_ROOT/memory", 'memory.limit_in_bytes' );
        }
        else {
            next;
        }

        my @parts = grep { $_ ne q{} } split m{/}xms, $path;
        for my $depth ( reverse 0 .. scalar @parts ) {
            my @ancestor = @parts[ 0 .. $depth - 1 ];
            push @files, File::Spec->catfile( $base, @ancestor, $file );
        }
    }

    return @files;
}

# Returns the limit in a cgroup memory limit file, in bytes, or nothing if the
# file does not exist or does not hold a number.
sub read_limit {
    my ($file) = @_;

    open my $fh, '<', $file
        or return;
    my $value = <$fh>;
    close $fh
        or die "Error: cannot close $file: $OS_ERROR\n";

    if ( !defined $value ) {
        return;
    }

    my ($limit) = $value =~ /\A (\d+) \s* \z/xms;
    if ( !defined $limit ) {
        return;
    }

    return $limit;
}

# Returns a field of /proc/meminfo in bytes, or undef if it is unavailable.
sub meminfo_bytes {
    my ($field) = @_;

    open my $fh, '<', $PROC_MEMINFO
        or return;
    my @lines = <$fh>;
    close $fh
        or die "Error: cannot close $PROC_MEMINFO: $OS_ERROR\n";

    for my $line (@lines) {
        my ($kib) = $line =~ /\A \Q$field\E : \s+ (\d+) \s+ kB/xms;
        if ( defined $kib ) {
            return $kib * $BYTES_PER_KIB;
        }
    }

    return;
}

# Calls STAR in order to build the index.
sub build_index {
    my ($arg) = @_;

    my $database        = $arg->{database};
    my $annotation_args = $arg->{annotation_args};
    my $memory          = $arg->{memory};

    my $name = basename( $database, qw(.fa .fas .fasta .fna) );
    report("STAR-indexing $name");

    make_path($INDEX_DIR);

    run('STAR',
        '--runThreadN'             => $THREADS,
        '--runMode'                => 'genomeGenerate',
        '--genomeDir'              => $INDEX_DIR,
        '--genomeFastaFiles'       => $database,
        '--limitGenomeGenerateRAM' => $memory,
        @{$annotation_args},
    );

    my $state
        = @{$annotation_args} ? 'with_annotation' : 'without_annotation';
    print {*STDOUT} "index $state\n"
        or die "Error: cannot write to standard output: $OS_ERROR\n";

    return;
}

# Stops the run if the index is too large for the memory STAR may use. Every
# alignment loads the whole index, so this is the least it will need.
sub check_index_fits {
    my ($star_memory) = @_;

    my $index_size
        = sum0 map { -s File::Spec->catfile( $INDEX_DIR, $_ ) || 0 }
        qw(Genome SA SAi);

    if ( $index_size > $star_memory ) {
        die 'Error: the genome index is '
            . gib($index_size)
            . ', which will not fit in the '
            . gib($star_memory)
            . ' that STAR may use; the job needs more than '
            . gib( $index_size + $MEMORY_HEADROOM )
            . " of memory\n";
    }

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

    # STAR rejects an option given twice, so a sort limit passed through by
    # the user takes precedence over the one the wrapper would add. STAR also
    # accepts the value after an equals sign or a space in the same argument.
    my @memory_args = ( '--limitBAMsortRAM' => $arg->{memory} );
    if ( any {/\A --limitBAMsortRAM (?: [=\s] | \z )/xms}
        @{ $arg->{star_args} } )
    {
        @memory_args = ();
    }

    my @align_command = (
        'STAR',
        @{ $arg->{annotation_args} },
        @{ $arg->{star_args} },
        '--runThreadN'        => $THREADS,
        '--genomeDir'         => $INDEX_DIR,
        '--outReadsUnmapped'  => 'Fastx',
        '--outFileNamePrefix' => "$prefix.",
        @memory_args,
        '--readFilesIn'      => @read_files,
        '--readFilesCommand' => 'gunzip',
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

# Formats a number of bytes as GiB.
sub gib {
    my ($bytes) = @_;

    return sprintf '%.2f GiB', $bytes / $BYTES_PER_GIB;
}

# Prints a message to stderr, exiting if the write fails.
sub report {
    my ($message) = @_;

    print {*STDERR} "$message\n"
        or die "Error: cannot write to standard error: $OS_ERROR\n";

    return;
}
