package orthomcl_core;
use strict;
use orthomcl_io;

require Exporter;

our @ISA= qw (Exporter);
our @EXPORT = qw(
              construct_graph
);

# this module contains 10 core subroutines used in orthomcl.pl
# they are:
# construct_graph
# make_inparalog
# make_ortholog
# rbh_weight
# satisfy_cutoff
# blastquery_ab
# non_redundant_list
# calc_weight
# calc_pmatch
# calc_matchlen


sub construct_graph {
    write_log("\nConstructing OrthoMCL similarity graph (for every pair of genomes)...\n");
    my %ortho;
    my $bpo_idx_ref;
    if ($::setting{'BPO_IDX_MODE'} eq 'all') {
        $bpo_idx_ref = bpo_idx_file('read','#all#');
    }
    foreach my $taxon (@::taxa) {
        write_log("\nIdentifying in-paralogs from $taxon\n");
        if ($::setting{'BPO_IDX_MODE'} eq 'taxon') {
            $bpo_idx_ref = bpo_idx_file('read',$taxon);
        }
        my ($edge_ref,$weight_ref,$sumw,$c) = make_inparalog($bpo_idx_ref,$taxon);                      # identification of inparalogs
        edge_file('write',$taxon,$taxon,$edge_ref);
        weight_file('write',$taxon,$taxon,$weight_ref);
    }
    for(my $i=0;$i<$#::taxa;$i++) {
        for(my $j=$i+1;$j<=$#::taxa;$j++) {
            write_log("\nIdentifying ortholog pairs between $::taxa[$i] and $::taxa[$j]\n");
            if ($::setting{'BPO_IDX_MODE'} eq 'taxon') {
                %{$bpo_idx_ref}=();
                my $bpo_idx_ref_i = bpo_idx_file('read',$::taxa[$i]);
                my $bpo_idx_ref_j = bpo_idx_file('read',$::taxa[$j]);
                foreach my $g (keys %$bpo_idx_ref_i) {
                    $bpo_idx_ref->{$g}=$bpo_idx_ref_i->{$g};
                }
                foreach my $g (keys %$bpo_idx_ref_j) {
                    $bpo_idx_ref->{$g}=$bpo_idx_ref_j->{$g};
                }
                undef $bpo_idx_ref_i; undef $bpo_idx_ref_j;
            }
            my ($edge_ref,$weight_ref,$sumw,$c_ortholog) = make_ortholog($bpo_idx_ref,$::taxa[$i],$::taxa[$j]); # identification of orthologs
            write_log("Appending co-ortholog pairs between $::taxa[$i] and $::taxa[$j]: ");  
            my $c_coortholog=0;
            my %edge_rbh;
            foreach my $pi (keys %$edge_ref) {@{$edge_rbh{$pi}}=@{$edge_ref->{$pi}};}  #make a copy of current %edge into %edge_rbh

            my %p1=%{edge_file('read',$::taxa[$i],$::taxa[$i])};
            my %p2=%{edge_file('read',$::taxa[$j],$::taxa[$j])};

            my %para;
            foreach my $p (keys %p1) {$para{$p}=$p1{$p};}
            foreach my $p (keys %p2) {$para{$p}=$p2{$p};}
            undef %p1; undef %p2;

            foreach my $n (keys %edge_rbh) {
                $ortho{$n} = 1;
                my (@nodes_1, @nodes_2);

                if (exists($para{$n})) {push (@nodes_1, $n, @{$para{$n}});}
                else {push (@nodes_1, $n);}

                foreach (@{$edge_rbh{$n}}) {
                    if (exists($para{$_})) {push (@nodes_2, $_, @{$para{$_}});}
                    else {push (@nodes_2, $_);}
                }

                @nodes_1=@{nonredundant_list(\@nodes_1)};
                @nodes_2=@{nonredundant_list(\@nodes_2)};

                foreach my $node_1 (@nodes_1) {
                    foreach my $node_2 (@nodes_2) {
                        next if(exists($weight_ref->{$node_1.' '.$node_2}));
                        my ($pv1,$flag1)=blastquery_ab($bpo_idx_ref,$node_1,$node_2);
                        my ($pv2,$flag2)=blastquery_ab($bpo_idx_ref,$node_2,$node_1);;
                        if (($flag1==1) && ($flag2==1)) {
                            push (@{$edge_ref->{$node_1}}, $node_2);
                            push (@{$edge_ref->{$node_2}}, $node_1);
                            my $w1=calc_weight($pv1);
                            my $w2=calc_weight($pv2);
                            my $w = ($w1+$w2)/2; # use averaged score as edge weight
                            $weight_ref->{$node_1.' '.$node_2} = sprintf("%.3f", $w);
                            $weight_ref->{$node_2.' '.$node_1} = sprintf("%.3f", $w);
                            $sumw += $w;
                            $c_coortholog++;
                        }
                    }
                }
            }
            write_log("$c_coortholog pairs\n");
            my $avgw = $sumw/($c_ortholog+$c_coortholog);
            write_log("$::taxa[$i] and $::taxa[$j] average weight: $avgw\n");
            foreach my $p (keys %$weight_ref) {
                $weight_ref->{$p} = sprintf("%.3f", $weight_ref->{$p}/$avgw);
            }

            edge_file('write',$::taxa[$i],$::taxa[$j],$edge_ref);
            weight_file('write',$::taxa[$i],$::taxa[$j],$weight_ref);
        }
    }

    undef $bpo_idx_ref;
    %::gindex=();

    foreach my $taxon (@::taxa) {
        write_log("\ncalculate average weight from $taxon\n");
        my %weight=%{weight_file('read',$taxon,$taxon)};

        my $count=0; my $sum=0;
        my $count_all=0; my $sum_all = 0;

        foreach my $pair (keys %weight) {
            my ($n,$p) = split(' ',$pair);
            $count_all++; $sum_all += $weight{$n.' '.$p};
            if ($ortho{$n} || $ortho{$p}) {
                $count++;
                $sum += $weight{$n.' '.$p};
            }
        }

        my $avgw = 0;
        
        # normalize the in-paralog weights by the average weight of inparalogs which have orthologs in other species
        # common case, for eukaryotes and most prokaryotes
        if ($count) {
            $avgw = $sum/$count;
        }
        # OR normalize the in-paralog weights by the average weight
        # not common
        elsif ($count_all) {
            $avgw = $sum_all/$count_all;
            write_log("taxon average weight is calculated based on all inparalog pairs\n");
        }
        # OR no normalization since $count_all=0 and there is nothing stored in %weight
        # not common, useful for prokaryotes or pathogens 

        write_log("$taxon average weight: $avgw\n");
        foreach my $p (keys %weight) {
            $weight{$p} = sprintf("%.3f", $weight{$p}/$avgw);
        }
        weight_file('write',$taxon,$taxon,\%weight);
    }
    %ortho=();

}

sub make_inparalog {
    my ($bpo_idx_ref,$taxon_id) = @_;
    my (%oneway_bh, %pvalue);
    foreach my $query_id (@{$::gindex{$taxon_id}}) {
        my ($offset_s,$offset_e);
        if (exists($bpo_idx_ref->{$query_id})) {
            ($offset_s,$offset_e)=split(";",$bpo_idx_ref->{$query_id});
        } else {next;}
        seek($::filehandle{'BPO'},$offset_s,0);
        my $pv_cut=10;
        LOOP:while (my $line=readline($::filehandle{'BPO'})) {
            $line=~s/\r|\n//g;
            my ($pv,$subject_id,$flag)=satisfy_cutoff($line);
            if ((tell($::filehandle{'BPO'})>$offset_e) || ($pv>$pv_cut)) {last LOOP;}
            next unless (($flag==1) && ($query_id ne $subject_id));
            if (not exists($::gindex2{$subject_id})) {
                write_log("$subject_id not defined in GG file:\n$line\n");
                next LOOP;
            }
            if ($::gindex2{$subject_id} ne $taxon_id) {$pv_cut=$pv; next LOOP;}
            push(@{$oneway_bh{$query_id}},$subject_id);
            $pvalue{$query_id.' '.$subject_id} = $pv;
        }
    }
    my $no_tmp = scalar(keys %oneway_bh);
    write_log("$no_tmp sequences have one-way better hits within species\n");
    return rbh_weight(\%oneway_bh, \%pvalue);
}

sub make_ortholog {
    my ($bpo_idx_ref,$ta,$tb) = @_;
    my (%oneway_bh,%pvalue);
    foreach my $query_id (@{$::gindex{$ta}}) {
        my ($offset_s,$offset_e);
        if (exists($bpo_idx_ref->{$query_id})) {
            ($offset_s,$offset_e)=split(";",$bpo_idx_ref->{$query_id});
        } else {next;}
        seek($::filehandle{'BPO'},$offset_s,0);
        my $pv_cut=10;
        LOOP:while (my $line=readline($::filehandle{'BPO'})) {
            $line=~s/\r|\n//g;
            my ($pv,$subject_id,$flag)=satisfy_cutoff($line);
            if ((tell($::filehandle{'BPO'})>$offset_e) || ($pv>$pv_cut)) {last LOOP;}
            if (not exists($::gindex2{$subject_id})) {
                write_log("$subject_id not defined in GG file:\n$line\n");
                next LOOP;
            }
            next unless (($flag==1) && ($::gindex2{$subject_id} eq $tb));
            push(@{$oneway_bh{$query_id}},$subject_id);
            $pvalue{$query_id.' '.$subject_id} = $pv;
            $pv_cut=$pv;
        }
    }
    my $no_tmpa=scalar(keys %oneway_bh);
    foreach my $query_id (@{$::gindex{$tb}}) {
        my ($offset_s,$offset_e);
        if (exists($bpo_idx_ref->{$query_id})) {
            ($offset_s,$offset_e)=split(";",$bpo_idx_ref->{$query_id});
        } else {next;}
        seek($::filehandle{'BPO'},$offset_s,0);
        my $pv_cut=10;
        LOOP:while (my $line=readline($::filehandle{'BPO'})) {
            $line=~s/\r|\n//g;
            my ($pv,$subject_id,$flag)=satisfy_cutoff($line);
            if ((tell($::filehandle{'BPO'})>$offset_e) || ($pv>$pv_cut)) {last LOOP;}
            if (not exists($::gindex2{$subject_id})) {
                write_log("$subject_id not defined in GG file:\n$line\n");
                next LOOP;
            }
            next unless (($flag==1) && ($::gindex2{$subject_id} eq $ta));
            push(@{$oneway_bh{$query_id}},$subject_id);
            $pvalue{$query_id.' '.$subject_id} = $pv;
            $pv_cut=$pv;
        }
    }
    my $no_tmpb=scalar(keys %oneway_bh)-$no_tmpa;

    write_log("$no_tmpa($ta)/$no_tmpb($tb) sequences have one-way best hits from the other species\n");
    return rbh_weight(\%oneway_bh, \%pvalue);
}

sub rbh_weight {
    my %oneway_bh = %{$_[0]};
    my %pvalue     = %{$_[1]};
    my (%edge, %weight);
    my $count=0;
    my $sumw=0;

    foreach my $query_id (sort keys %oneway_bh) {
        foreach my $subject_id (@{$oneway_bh{$query_id}}) { # all the subject_id is oneway best hit of query_id
            next if($weight{$query_id.' '.$subject_id}); # this pair was already identified as RBH pair

            if (exists($pvalue{$subject_id.' '.$query_id})) { # means query_id is also oneway best hit of subject_id, so RBH pair
                push (@{$edge{$query_id}}, $subject_id);
                push (@{$edge{$subject_id}}, $query_id);
                write_rbh("$query_id    $subject_id    ".$pvalue{$query_id.' '.$subject_id}."    ".$pvalue{$subject_id.' '.$query_id}."\n");
                my $w1=calc_weight($pvalue{$query_id.' '.$subject_id});
                my $w2=calc_weight($pvalue{$subject_id.' '.$query_id});
                my $w = ($w1+$w2)/2;
                $sumw += $w;
                $count++;
                # use averaged score as edge weight
                $weight{$query_id.' '.$subject_id} = sprintf("%.3f", $w);
                $weight{$subject_id.' '.$query_id} = sprintf("%.3f", $w);
            }
        }
    }
    my $no_tmp = scalar(keys %weight)/2;
    write_log("$no_tmp sequence pairs were identified as Reciprocal Better/Best Hit\n");
    return (\%edge, \%weight, $sumw, $count);
}

sub calc_weight {
    #use -log10(P) as weights and treat P=0 as -log10(P)=$::setting{'MAX_WEIGHT'}
    my $pvalue=$_[0];
    if($pvalue == 0) {
        return $::setting{'MAX_WEIGHT'};
    } else {
        return -log($pvalue)/log(10);
    }
}

sub satisfy_cutoff {
#1;At1g01190;535;At1g01190;535;0.0;97;1:1-535:1-535.
#0 1         2   3         4   5   6  7
    my ($subject_id,$pv,$pi,$hsp_info,$query_len,$subject_len)=(split(";",$_[0]))[3,5,6,7,2,4];
    my $flag=1;   # 1, satisfy cutoff; 0, otherwise.
    
    if (exists($::setting{'PVALUE_CUTOFF'})) {
        if($pv > $::setting{'PVALUE_CUTOFF'}) {
            $flag=0;
            return ($pv,$subject_id,$flag);
        }
    }
    if (exists($::setting{'PIDENT_CUTOFF'})) {
        if ($pi < $::setting{'PIDENT_CUTOFF'}) {
            $flag=0;
            return ($pv,$subject_id,$flag);
        }
    }
    if (exists($::setting{'PMATCH_CUTOFF'})) {
        if (calc_pmatch($query_len,$subject_len,$hsp_info) < $::setting{'PMATCH_CUTOFF'}) {
            $flag=0;
            return ($pv,$subject_id,$flag);
        }
    }
    return ($pv,$subject_id,$flag);
}


sub blastquery_ab {
    my ($bpo_idx_ref,$a,$b) = @_;

    my ($offset_s,$offset_e);
    if (exists($bpo_idx_ref->{$a})) {
        ($offset_s,$offset_e)=split(";",$bpo_idx_ref->{$a});
    } else {return (undef,0);}
    seek($::filehandle{'BPO'},$offset_s,0);
    while (my $line=readline($::filehandle{'BPO'})) {
        $line=~s/\r|\n//g;
        my ($pv,$subject_id,$flag)=satisfy_cutoff($line);
        if (tell($::filehandle{'BPO'})>$offset_e) {return (undef,0);}
        next unless (($flag==1) && ($subject_id eq $b));
        return ($pv,$flag);
    }
    return (undef,0);
}

sub nonredundant_list {
    my $list_ref=$_[0];
    my %nr;
    foreach (@{$list_ref}) {$nr{$_}=1;}
    my @nr=sort (keys %nr);
    return \@nr;
}

sub calc_pmatch {
    my ($q_len,$s_len,$hsp_info) = @_;
    my (%s_start, %s_length, %q_start, %q_length);
    my @hsp=split(/\./,$hsp_info);
    foreach (@hsp) {
        if (/(\d+)\:(\d+)\-(\d+)\:(\d+)\-(\d+)/) {
            $s_start{$1}=$4; 
            $s_length{$1}=$5-$4+1;
            $q_start{$1}=$2;
            $q_length{$1}=$3-$2+1;
        }
    }
    my $s_matchlen = calc_matchlen(\%s_start,\%s_length);
    my $q_matchlen = calc_matchlen(\%q_start,\%q_length);
    if ($s_len >= $q_len) {
        return 100*$q_matchlen/$q_len;
    }else{
        return 100*$s_matchlen/$s_len;
    }
}

sub calc_matchlen {
    my %start        = %{$_[0]}; 
    my %length       = %{$_[1]};
    my @starts = sort{$start{$a}<=>$start{$b}} (keys %start);
    return $length{$starts[0]} if(scalar(@starts)==1);
    my $i=1; 
    my $pos =  $start{$starts[0]} + $length{$starts[0]};
    my $match_length = $length{$starts[0]}; 
    while ($i<scalar(@starts)) {
        if ($length{$starts[$i]} + $start{$starts[$i]} <= $pos) {
            $i++;
            next;
        }
        if ($start{$starts[$i]} > $pos) {
            $match_length += $length{$starts[$i]};
            $pos = $start{$starts[$i]} + $length{$starts[$i]};
        } else {
            $match_length += $length{$starts[$i]} - ($pos - $start{$starts[$i]});
            $pos = $start{$starts[$i]} + $length{$starts[$i]};
        }
        $i++;
    }

    return $match_length;
}


1;