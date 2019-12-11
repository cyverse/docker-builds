#!/bin/bash

export PERL5LIB=/etc/perl
export PATH=$PATH:/usr/bin

_usage() {
echo "##### U S A G E : Help and ERROR ######"
echo "$0 $Options"
echo "*"
echo "       Usage: $0 <[options]>"
echo "       Options:"
echo "               -h   --help             Show this message"
echo "               --gg_file               Path to input GG file created by fastaRename app."
echo "               --bpo_file              Path to parsed BLAST output created by parseBlastBpo app." 
echo "               --bpo_idx_mode          Advanced Use Only.  Mode to index BPO file. Options: 'all' (DEFAULT) or 'XX' taxon abbreviation from GG file."
echo "               --pvalue_cutoff         Advanced Use Only.  Equivalent to BLAST e-value cutoff. BLAST similarities with an e-value greater than this value are ignored during clustering.  1e-5 (DEFAULT)."
echo "               --pident_cutoff         Advanced Use Only.  Equivalent to BLAST percent identity cutoff. BLAST similarities with percent identity lower than this value are ignored during clustering.  <0-100>, 0 (DEFAULT)"
echo "               --pmatch_cutoff         Advanced Use Only.  Percent match cutoff <0-100> in ortholog clustering, 0 (DEFAULT)"
echo "               --max_weight            Advanced Use Only.  Weight used for protein pairs whose BLASTp p-value is zero. This depends on the algorithm you use: if the second smallest p-value is in the order of -99, maximum_weight should be 100; if -299, maximum_weight should be 300 <DEFAULT>"
echo "               --inflation             Advanced Use Only.  Markov Inflation Index, used in MCL algorithm. Increasing this index increases cluster tightness, and the number of clusters."
echo "               --project_run_name      Will be the name of tne conf and log files"
echo ""
if [ "$1" != "" ]; then
echo "USAGE ERROR: $1"
fi
exit
}


main() {

    Options=$@
    OptNum=$#

    ## Initialize command line vars
    #branch_name=""
    #working_dir=""
    #virtualenv_dir=""
    #ssh_key_dir=""
    #setup_files_dir=""
    #server_name=""

    gg_file=""
    bpo_file=""
    bpo_idx_mode=""
    pvalue_cutoff=""
    pident_cutoff=""
    pmatch_cutoff=""
    max_weight=""
    inflation=""
    project_run_name=""
     
    ## Collect command line vars
    while getopts 'ahtT-b:D:E:s:' OPTION ; do
      case "$OPTION" in
        h  ) _usage                    ;;
        -  ) [ $OPTIND -ge 1 ] && optind=$(expr $OPTIND - 1 ) || optind=$OPTIND
             eval OPTION="\$$optind"
	     #echo "New arguement"
             #echo "OPTION: $OPTION"
	     OPTARG=$(echo $OPTION | cut -d'=' -f2)
             #echo  "OPTARG: $OPTARG"
	     OPTION=$(echo $OPTION | cut -d'=' -f1)
	     #echo "OPTION: $OPTION"
             case $OPTION in
                 --gg_file ) gg_file="$OPTARG"             ;;
                 --bpo_file ) bpo_file="$OPTARG"           ;;
                 --bpo_idx_mode) bpo_idx_mode="$OPTARG"    ;;
                 --pvalue_cutoff ) pvalue_cutoff="$OPTARG" ;;
                 --pident_cutoff ) pident_cutoff="$OPTARG" ;;
                 --pmatch_cutoff ) pmatch_cutoff="$OPTARG" ;;
                 --max_weight ) max_weight="$OPTARG"       ;;
                 --help ) _usage                           ;;
                 --inflation ) inflation="$OPTARG"         ;;
		 --project_run_name ) project_run_name="$OPTARG"  ;;
                 * )  _usage " Unknown long option provided: $OPTION:$OPTARG" ;;
             esac
           OPTIND=1
           shift
          ;;
        ? )  _usage " Unknown short option provided: $OPTION:$OPTARG" ;;
      esac
    done


   
}

build_conf_file(){
  compare_string=""


  filename="$project_run_name.conf"
  echo "filename: $filename"
  echo "project: $project_run_name"
  current_working_directory=`pwd`;
  touch "$filename"
  echo "DIR=$current_working_directory/" >> $filename
  echo "GG_FILE=$gg_file" >> $filename
  echo "BPO_FILE=$bpo_file" >> $filename
  echo "BPO_IDX_MODE=$bpo_idx_mode" >> $filename

  if [ "$pvalue_cutoff" != "$compare_string" ]; then
     echo "PVALUE_CUTOFF=$pvalue_cutoff" >> $filename
  fi

  if [ "$pident_cutoff" != "$compare_string" ]; then
     echo "PIDENT_CUTOFF=$pident_cutoff" >> $filename
  fi

  if [ "$pmatch_cutoff" != "$compare_string" ]; then
    echo "PMATCH_CUTOFF=$pmatch_cutoff" >> $filename
  fi
 
  echo "MAX_WEIGHT=$max_weight" >> $filename
  echo "MCL=/usr/local/bin/mcl" >> $filename
  echo "INFLATION=$inflation" >> $filename
}

run(){
  touch "$project_run_name.log"
  /usr/bin/orthomcl.pl "$project_run_name.conf" > "$project_run_name.log" 2>&1
}

#EXECUTION PATH:
main "$@"
build_conf_file
run
