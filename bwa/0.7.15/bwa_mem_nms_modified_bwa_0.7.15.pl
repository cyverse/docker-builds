#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

#my $app = "bwa";
my $app = "/bwa-0.7.15/bwa";

# Define worflow options
my ($query_file1, $query_file2, $database_path, $user_database_path);

my $null;
my $format = 'SE';

GetOptions( "query|query1|input=s" => \$query_file1,
                       "query2=s" => \$query_file2,
                       "database=s" => \$database_path,
                       "user_database=s" => \$user_database_path,
                       "output=s" => \$null );

# Allow over-ride of system-level database path with user
# May not need to do this going forward...
if (defined($user_database_path)) {
       $database_path = $user_database_path;
}

my $out_file = "bwa_output.sam";

# Grab any flags or options we don't recognize and pass them as plain text
# Need to filter out options that are handled by the GetOptions call
my @args_to_reject = qw(-b -c);
my $BWA_ARGS = join(" ", @ARGV);
foreach my $a (@args_to_reject) {
       if ($BWA_ARGS =~ /$a/) {
               report("Most bwa arguments are legal for use with this script, but $a is not. Please omit it and submit again");
               exit 1;
       }
}

# Check for presence of second read file
if (defined($query_file2)) {
       $format = 'PE';
       $out_file = $query_file1 . "-" . $query_file2 . ".sam";
       report("Pair-end alignment requested");
}

# Auto index
unless(-e "$database_path.bwt"){
       my $fsize = -s $database_path;
       my $algo = "is";
       if($fsize > 100000000){
           $algo = "bwtsw";
       }
       report("Indexing using $algo");
       system("$app index -a $algo $database_path");
}

unless (-e "$database_path.bwt") {
       report("Indexing failure...");
       exit 1;
}

my @stack;
#push(@stack, "$app mem $BWA_ARGS $database_path $query_file1 > temp1.sai");

if ($format eq 'SE') {
#	push(@stack, "$app samse $database_path temp1.sai $query_file1 > $out_file");
	push(@stack, "$app mem $BWA_ARGS $database_path $query_file1 > $out_file");
} elsif ($format eq 'PE') {
#       push(@stack, "$app aln $BWA_ARGS $database_path $query_file2 > temp2.sai");
        push(@stack, "$app mem $BWA_ARGS $database_path $query_file1 $query_file2 > $out_file");
}

foreach my $cmd (@stack) {
       report($cmd);
       system($cmd);
}

if (-e $out_file) {
       exit 0
} else {
       exit 1
}


sub report {

       my $text = shift;
       print STDERR $text, "\n";

}
