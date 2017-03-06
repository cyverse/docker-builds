#!/usr/bin/perl -w
use strict;
use File::Copy qw/move copy/;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper;

use constant CUFFLINKS  => ('2.2.1' => '/cufflinks-2.2.1.Linux_x86_64/');

report_input_stack();

# Define workflow options
my (@query_file, @query_file1, $annotation, $user_annotation, $version, $contained);
$version = '2.2.1';

GetOptions( "infile=s"    => \@query_file,
	    "G=s"         => \$annotation,
            "M=s"         => \$user_annotation,
	    "version=s"   => \$version,
	    "contained"   => \$contained,
	    );

# more with this later, when they fix BOOLEANs

my (@queries,$success);

if (@query_file) {
    push @queries, @query_file;
}


@queries > 0 || die "I could not find any GTF input files.\n";


# Allow over-ride of system-level database path with user
# May not need to do this going forward...
if (defined($user_annotation)) {
    $annotation = $user_annotation;
    `dos2unix $user_annotation`;
}

# Grab any flags or options we don't recognize and pass them as plain text
# Need to filter out options that are handled by the GetOptions call
my @args_to_reject = qw(-xxxx);
my $CUFFCOMPARE_ARGS = join(" ", @ARGV);
foreach my $a (@args_to_reject) {
    if ($CUFFCOMPARE_ARGS =~ /$a/) {
	report("Most arguments are legal for use with this script, but $a is not. Please omit it and submit again");
	exit 1;
    }
}


my %app = CUFFLINKS;
my $clpath = $app{$version} or die "Cufflinks version $version is not supported."; 
my $cuffcompare = $clpath . "cuffcompare -o cuffcompare_out";
my $cmd = $annotation ? "$cuffcompare $CUFFCOMPARE_ARGS $annotation " : "$cuffcompare $CUFFCOMPARE_ARGS";
chomp($ENV{PATH} = `echo \$PATH`);
$ENV{PATH} = join(':',$ENV{PATH},$clpath);




#my $gtf_out = 'gtf';
for my $query_file (@queries) {    
    my $basename = $query_file;
    $basename =~ s/^\S+\/|\.\S+$//g;

    my $cuffcommand = $cmd . " $query_file";
    report("Executing: $cuffcommand");

    system("$cuffcommand");
    $success++ if -e "cuffcompare_out.combined.gtf" && ! -z "cuffcompare_out.combined.gtf";
}


die "Something did not work, no combined.gtf file!" unless $success;

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

