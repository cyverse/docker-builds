#!/usr/bin/perl -w
use strict;


while (<>) {
    chomp;
    my ($class_code) = /class_code "([=cj])"/;
    if ($class_code) {
        my ($tid) = /oId "([^\"]+)"/;
        if ($tid) {
            (my $gid = $tid) =~ s/\.\S+$//;
            s/gene_id "[^\"]+"/gene_id "$gid"/;
            s/transcript_id "[^\"]+"/transcript_id "$tid"/;
	    s/; oId.+class_code "=";//;
            print "$_\n";
        }
        else {
            print "$_\n";
        }
    }
    else {
        print "$_\n";
    }
}

