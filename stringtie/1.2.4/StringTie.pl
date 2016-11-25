#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


# remember to remove this
report_input_stack();

my (@file_query, $user_annotation_path, $annotation_path);

GetOptions( "file_query=s"      => \@file_query,
            "user_annotation=s" => \$user_annotation_path
	    );

if ($user_annotation_path) {
    $annotation_path = $user_annotation_path;
}


my $success = undef;



for my $query_file (@file_query) {
    # Grab any flags or options we don't recognize and pass them as plain text
    # Need to filter out options that are handled by the GetOptions call
    my @args_to_reject = qw(-xxxx);

    my $STRINGTIE_ARGS = join(" ", @ARGV);
    foreach my $a (@args_to_reject) {
	if ($STRINGTIE_ARGS =~ /$a/) {
	    report("Most TopHat arguments are legal for use with this script, but $a is not. Please omit it and submit again");
	    exit 1;
	}
    }
   my $stringtie = "stringtie"; 
	$STRINGTIE_ARGS .= " -G $annotation_path";
  chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;
    system "mkdir $basename";
    my $align_command = "$stringtie $query_file -o $basename.gtf -B $STRINGTIE_ARGS -C $basename.refs.gtf -A $basename.abund.tab";
    
    report("Executing: $align_command\n");
    system $align_command;
    system "mv $basename.gtf $basename.refs.gtf *.abund.tab *ctab $basename";
}

sub report {
    print STDERR "$_[0]\n";
}

sub report_input_stack {
    my @stack = @ARGV;
    report(Dumper \@stack);
}
