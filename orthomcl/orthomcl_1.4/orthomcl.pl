#!/usr/bin/perl -w
use strict;
use orthomcl_io;
use orthomcl_core;


our (%setting,%filehandle);

our (@taxa,%gindex,%gindex2,@mcl_idx);
our (%edge_all,%weight_all);

if (!defined $ARGV[0]) {die "SETTING_FILE needs to be given?";}

$setting{'START_TIME'}=`date`;
read_setting($ARGV[0]);
%filehandle=%{setup_run()};
write_log("OrthoMCL Clustering Starts:\nRUN DIRECTORY: $setting{'RUN_DIR'}\n");
read_ggfile();
if (exists($setting{'FORMER_RUN_DIR'})) {
    @taxa=qw($setting{'TAXA_LIST'}) if (exists($setting{'TAXA_LIST'}));
    write_log("Constructing the whole graph (reading from $setting{FORMER_RUN_DIR})...\n");
} else {
    index_bpo();
    construct_graph();
    write_log("Constructing the whole graph...\n");
}

for(my $i=0;$i<=$#taxa;$i++) {
    for(my $j=$i;$j<=$#taxa;$j++) {
        my %edge   = %{edge_file('read',$taxa[$i],$taxa[$j])};
        my %weight = %{weight_file('read',$taxa[$i],$taxa[$j])};
        foreach my $n (keys %edge) {
            push(@{$edge_all{$n}},@{$edge{$n}});
        }
        foreach my $p (keys %weight) {
            $weight_all{$p} = $weight{$p};
        }
    }
}

write_matrix_index();
%edge_all=();
%weight_all=();
execute_MCL();
mcl_backindex();

$setting{'END_TIME'}=`date`;
write_log("\nStart Time: $setting{'START_TIME'}\nEnd Time:   $setting{'END_TIME'}\n");
write_setting();

