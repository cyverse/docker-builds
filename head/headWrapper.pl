#!/usr/bin/perl

use strict;
use warnings;

# a wrapper for head for use in galaxy
# headWrapper.pl [filename] [# lines to show] [output]

die "Check arguments" unless @ARGV == 3;

open (OUT, ">$ARGV[2]") or die "Cannot create $ARGV[2]:$!\n";
open (HEAD, "head -n $ARGV[1] $ARGV[0]|") or die "Cannot run head:$!\n";
while (<HEAD>) {
    print OUT;
}
close OUT;
close HEAD;
