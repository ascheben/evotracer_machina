#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage: run_pipeline.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--infile) ASV="$2"; shift ;;
        -t|--tree) TREE="$2"; shift ;;
        -s|--scripts) SPATH="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; echo "Usage: run_pipeline.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts>" ; exit 1 ;;
    esac
    shift
done

if ! command -v pmh_tr &> /dev/null
then
    echo "pmh_tr from the MACHINA installation could not be found. Exiting!"
    exit
fi

# PATHS TO SCRIPTS 

TRAV="${SPATH}/traverse_split.py"
GET="${SPATH}/get_results.sh"
MACHINA="pmh_tr"
TOPOLOGY="${SPATH}/print_seeding_topology.py"
SELECTION="${SPATH}/selection_test.py"
MIGRATION="${SPATH}/count_migrations.py"

## PREPROCESS INPUT DATA ##

# Extract key ASV columns
#asv_names,sample,group
cut -d',' -f1,2,30 ${ASV} > asv_sample_group.csv
# exclude miscelleneaous group CP00
cut -f3 -d',' asv_sample_group.csv|tail -n +2| sort| uniq| grep -v CP00| grep -v CP01| grep -v CP02| grep -v CP03 > CP_list.txt
# make output directory for each CP
while read l; do mkdir ${l}_split;done<CP_list.txt
# prepare machina input files for each CP
while read l; do ${TRAV} tree_all_clones.newick asv_sample_group.csv ${l};done<CP_list.txt

## RUN MACHINA ##     

# make command file for machina
while read l; do echo "${MACHINA} -o ${l}_split -c ${l}_colors.txt -p PRL ${l}_tree_split.txt ${l}_labels_split.txt &> ${l}_split/results.txt";done<CP_list.txt > machina.cmd
# run machina in parallel
ParaFly -CPU 12 -c machina.cmd
# parse results from each machina output dir
while read l; do ${GET} $l;done<CP_list.txt | tr '\t' ' '>> all_results.txt
# move intermediate output
mkdir cp_output
mv CP* cp_output
mv machina* cp_output
mv asv_sample_group.csv cp_output

## ANALYSE INFERRED TOPOLOGY



python $TOPOLOGY all_results.txt > seeding_topology.txt 
python $SELECTION all_results.txt asv_stat.csv > selection.txt
python $MIGRATION all_results.txt > migration.txt

