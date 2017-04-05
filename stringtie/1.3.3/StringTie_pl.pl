#!/usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);


# remember to remove this
report_input_stack();

my (@file_query,$STRINGTIE_ARGS);

GetOptions( "file_query=s"      => \@file_query
	    );
system "mkdir StringTie_output gtf_files ballgown_input_files";
$STRINGTIE_ARGS = join(" ", @ARGV);
for my $query_file (@file_query) {
    
   my $stringtie = "stringtie"; 
  chomp(my $basename = `basename $query_file`);
    $basename =~ s/\.\S+$//;
    print "$STRINGTIE_ARGS\n";
    if ($STRINGTIE_ARGS =~ /-G/) {
    my $align_command = "$stringtie $query_file -o $basename.gtf $STRINGTIE_ARGS -C $basename.refs.gtf -A $basename.abund.tab";
    print "using reference annotation\n";
    report("Executing: $align_command\n");
    system $align_command;
    }
    else{
        my $align_command = "$stringtie $query_file -o $basename.gtf $STRINGTIE_ARGS ";
        print "no annotation used\n";
    report("Executing: $align_command\n");
    system $align_command;
    }
    system "mkdir StringTie_output/$basename ballgown_input_files/$basename";
    system "mv $basename.gtf $basename.refs.gtf *.abund.tab *ctab StringTie_output/$basename";
    system "rsync -av --progress StringTie_output/$basename/$basename.gtf gtf_files";
    system "rsync -av --progress StringTie_output/$basename/$basename.gtf StringTie_output/$basename/*ctab ballgown_input_files/$basename";
}

sub report {
    print STDERR "$_[0]\n";
}

sub report_input_stack {
    my @stack = @ARGV;
    report(Dumper \@stack);
}
