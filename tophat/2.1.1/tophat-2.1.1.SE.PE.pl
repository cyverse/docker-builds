#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use constant TOPHATP    => ('2.1.1'  => '/tophat-2.1.1.Linux_x86_64/');
use constant BOWTIEP   => ('2.2.5'  => '/bowtie2-2.2.5/');

use constant SAMTOOLSP  => '/usr/bin/samtools/';

# remember to remove this
report_input_stack();

my (@file_query, $database_path, $user_database_path, $annotation_path, 
$user_annotation_path, $file_names, $root_names, @file_query2,$together,$colorspace);

my $version   = '2.1.1';
my $btversion = '2.2.5';

GetOptions( "file_query=s"      => \@file_query,
	    "file_query2=s"     => \@file_query2,
	    "database=s"        => \$database_path,
	    "user_database=s"   => \$user_database_path,
	    "annotation=s"      => \$annotation_path,
	    "user_annotation=s" => \$user_annotation_path,
	    "file_names=s"      => \$file_names,
	    "root_names=s"      => \$root_names,
	    "tophat_version=s"  => \$version,
	    "bowtie_version=s"  => \$btversion,
	    "align_reads=s"     => \$together,
	    "usecolorspace"	=> \$colorspace
	    );

# Enabling colorspace support
# Die if bowtie version is > 2 and we have specified colorspace
if (($colorspace) and ($btversion =~ /^2/)) {
	die "You have requested to align in colorspace but asked to use Bowtie 2.x. Please resubmit requesting usage of Bowtie 0.12.7";
}

# sanity check for input data
my $pe;
if (@file_query2) {
    @file_query && @file_query2 || die "Error: At least one file for each paired-end is required\n"; 
    @file_query == @file_query2 || die "Error: Unequal number of files for paired ends\n";
    $pe++;
}

unless ($together eq 'together') {
    undef $together;
}


if (!($user_database_path || $database_path)) {
    die "No reference genome was supplied\n";
}
if (@file_query < 1) {
    die "No FASTQ files were supplied\n";
}

my %tophat = TOPHATP;
my %bowtie = BOWTIEP;

my $tophatp   = $tophat{$version};
my $bowtiep   = $bowtie{$btversion};

chomp($ENV{PATH} = `echo \$PATH`);
$ENV{PATH} = join(':',$ENV{PATH},$tophatp,$bowtiep,SAMTOOLSP);


# Sanity check for input ref. genome
unless ($database_path || $user_database_path) {
  die "No reference genome was selected" 
}

# Allow over-ride of system-level database path with user
if ($user_database_path) {
  $database_path = $user_database_path;
  unless (`grep \\> $database_path`) {
      die "Error: $database_path  the user supplied file is not a FASTA file";
  }
  my $name = basename($database_path, qw/.fsa .fa .fas .fasta .fna/);
  print STDERR "bowtie-indexing $name\n";
  $bowtiep .= $btversion =~ /^2/ ? 'bowtie2-build' : 'bowtie-build';

  # rename genome file to use the .fa extension
  if ($database_path !~ /$name\.fa$/) {
      my $new_path = $database_path;
      $new_path =~ s/$name\.\S+$/$name\.fa/;
      system "cp $database_path $new_path";
  }
  
  # Enable bowtie-build to index in colorspace
  my $btidxarg = "";
  if ($colorspace) {
  	$btidxarg = "--color "
  }
  system $bowtiep . "$btidxarg $database_path $name";
  $database_path = $name;
}
if ($user_annotation_path) {
    $annotation_path = $user_annotation_path;
}

# Should be a temporary hack
if ($btversion =~ /^2/ && !$user_database_path) {
    my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
    unless (-e "$database_path.1.bt2") {  
	print STDERR "bowtie-indexing $name\n";
	system $bowtiep . "bowtie2-build $database_path $name";
	$database_path = $name;
    }
}

# If we are using new tophat with old bowtie
if ($version =~ /^2/ && $btversion eq '0.12.7') {
    push @ARGV, '--bowtie1';
}

my $success = undef;

my (@basenames,%sample);


my $nocount;
my $samples;
my @to_move;# = ('bam');

if ($together) {
    my $file_string1 = join(',',@file_query);
    @file_query = ($file_string1);
    if ($pe) {
	my $file_string2 = join(',',@file_query2);
	@file_query2 = ($file_string2);
    }
}

for my $query_file (@file_query) {
    # Grab any flags or options we don't recognize and pass them as plain text
    # Need to filter out options that are handled by the GetOptions call
    my @args_to_reject = qw(-xxxx);

    my $second_file;
    if ($pe) {
	$second_file = shift @file_query2 || die "Right reads file is missing!\n";
    }

    my $TOPHAT_ARGS = join(" ", @ARGV);
    foreach my $a (@args_to_reject) {
	if ($TOPHAT_ARGS =~ /$a/) {
	    report("Most TopHat arguments are legal for use with this script, but $a is not. Please omit it and submit again");
	    exit 1;
	}
    }

    my $app  = $tophatp.'tophat';
    if ($annotation_path) {
		$TOPHAT_ARGS .= " -G $annotation_path";
    }
    
    # Enable colorspace TopHat query
    if ($colorspace) {
    	$TOPHAT_ARGS .= " --color"
    }
    
    $second_file ||= '';
    my $align_command = "$app $TOPHAT_ARGS $database_path $query_file $second_file ";
    
    chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;

    report("Executing: $align_command\n");
    system $align_command;
    die "\n\nTopHat run failed: no accepted_hits.bam file\n" unless  -e "tophat_out/accepted_hits.bam";

    #$success or next;
    unless ($together) {
    	system "mv tophat_out $basename\_out";
	my $bam = 'bam';
	mkdir($bam) unless -d $bam;
	system "cp $basename\_out/accepted_hits.bam bam/$basename.bam";
	system "samtools index bam/$basename.bam";
    }
}
system "rm -f *.bt2";

#$success ? exit 0 : exit 1;

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
