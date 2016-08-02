#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


#use constant Sailfish    => ('Beta-0.9.0'  => '~/SailfishBeta-0.9.0/bin/');

# remember to remove this
#report_input_stack();

my (@file_query, $database_path, $user_database_path, $annotation_path, 
$user_annotation_path, $file_names, $root_names, @file_query2, $file_type, $lib_type);


GetOptions( "file_query=s"      => \@file_query,
	    "file_query2=s"     => \@file_query2,
	    "database=s"        => \$database_path,
	    "user_database=s"   => \$user_database_path,
            "annotation=s"      => \$annotation_path,
            "user_annotation=s" => \$user_annotation_path,
	    "file_names=s"      => \$file_names,
	    "root_names=s"      => \$root_names,
	    "file_type=s"       => \$file_type,
	    "lib_type=s"        => \$lib_type,
	    );

# sanity check for input data
if (@file_query2) {
    @file_query && @file_query2 || die "Error: At least one file for each paired-end is required\n"; 
    @file_query == @file_query2 || die "Error: Unequal number of files for paired ends\n";
}

if (!($user_database_path || $database_path)) {
    die "No reference set of transcripts was supplied\n";
}
if (@file_query < 1) {
    die "No FASTQ files were supplied\n";
}


# Allow over-ride of system-level database path with user
my $sailfish  = "/sailfish/bin/sailfish";

if ($user_database_path) {
  $database_path = $user_database_path;
  unless (`grep \\> $database_path`) {
      die "Error: $database_path the user supplied file is not a FASTA file";
  }
  my $name = basename($database_path, qw/.fa .fas .fasta .fna/);
  print STDERR "sailfish-indexing $name\n";
  system $sailfish . " index -t $database_path -o index";
  
  if ($database_path !~ /$name\.fa$/) {
      my $new_path = $database_path;
      $new_path =~ s/$name\.\S+$/$name\.fa/;
      #system "cp $database_path $new_path";
  }
  $database_path = $name;
}


my $success = undef;

#system "mkdir output";

for my $query_file (@file_query) {
    # Grab any flags or options we don't recognize and pass them as plain text
    # Need to filter out options that are handled by the GetOptions call
    my @args_to_reject = qw(-xxxx);


    my $second_file = shift @file_query2 if @file_query2;

    my $SAILFISH_ARGS = join(" ", @ARGV);
    foreach my $a (@args_to_reject) {
	if ($SAILFISH_ARGS =~ /$a/) {
	    report("Most Sailfish arguments are legal for use with this script, but $a is not. Please omit it and submit again");
	    exit 1;
	}
    }

# Check for presence of second read file
#if (defined($second_file)) {
#       $format = 'PE';
#       report("Pair-end alignment requested");
#}

my $format = $file_type;
my $library = $lib_type;


chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;

if ($format eq 'PE') {
          my $align_command = "$sailfish quant -i index -l $library -1 $query_file -2 $second_file -o $basename $SAILFISH_ARGS";
          system $align_command;
                     } 
elsif($format eq 'SE'){
          my $align_command = "$sailfish quant -i index -l $library -r $query_file -o $basename $SAILFISH_ARGS";
          system $align_command;        
        	}
   	}


sub report {
    print STDERR "$_[0]\n";
}

sub report_input_stack {
    my @stack = @ARGV;
    report(Dumper \@stack);
}
