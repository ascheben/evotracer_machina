#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo "Usage: run_pipeline.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts> --prefix <output_prefix>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--infile) ASV="$2"; shift ;;
        -t|--tree) TREE="$2"; shift ;;
        -s|--scripts) SPATH="$2"; shift ;;
        -p|--prefix) PREFIX="$2"; shift ;;
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
SELECTION="${SPATH}/selection_tree_test.py"
MIGRATION="${SPATH}/count_migrations.py"

## PREPROCESS INPUT DATA ##

# Extract key ASV columns
#asv_names,sample,group
cut -d',' -f1,2,30 ${ASV} > ${PREFIX}_asv_sample_group.csv
# exclude miscelleneaous group CP00
cut -f3 -d',' ${PREFIX}_asv_sample_group.csv|tail -n +2| sort| uniq| grep -v CP00 > ${PREFIX}_CP_list.txt
# prepare machina input files for each CP
while read l; do ${TRAV} ${TREE} ${PREFIX}_asv_sample_group.csv ${l};done<${PREFIX}_CP_list.txt
# remove CPs that failed traversal due to unexpected clade topology
#grep -v -f FailedCP.txt ${PREFIX}_CP_list.txt > ${PREFIX}_CP_list.tmp
#mv ${PREFIX}_CP_list.tmp ${PREFIX}_CP_list.txt
# make output directory for each CP
while read l; do mkdir ${l}_split;done<${PREFIX}_CP_list.txt
   
# prepare input for selection scan on original tree
mkdir ${PREFIX}_cp_output
for l in *_tree_split.txt; do cat $l |sed "s/^/${l} tree /"| sed 's/_tree_split.txt//';done | grep -v CP00 > ${PREFIX}_cp_output/${PREFIX}_all_original_tree.txt

## RUN MACHINA ##     
# make command file for machina
# speed up MACHINA a bit by adding "-m 3"
while read l; do echo "${MACHINA} -m 1 -o ${l}_split -c ${l}_colors.txt -p PRL ${l}_tree_split.txt ${l}_labels_split.txt &> ${l}_split/results.txt";done<${PREFIX}_CP_list.txt > ${PREFIX}_machina.cmd
# run machina in parallel
ParaFly -CPU 12 -c ${PREFIX}_machina.cmd
# parse results from each machina output dir
while read l; do ${GET} $l;done<${PREFIX}_CP_list.txt | tr '\t' ' '>> ${PREFIX}_cp_output/${PREFIX}_all_results.txt
# move intermediate output
mv CP* ${PREFIX}_cp_output
mv ${PREFIX}_machina* ${PREFIX}_cp_output
mv ${PREFIX}_asv_sample_group.csv ${PREFIX}_cp_output
#mv FailedCP.txt ${PREFIX}_cp_output
mv ${PREFIX}_CP_list.txt ${PREFIX}_cp_output
## ANALYSE INFERRED TOPOLOGY

python $TOPOLOGY ${PREFIX}_cp_output/${PREFIX}_all_results.txt > ${PREFIX}_cp_output/${PREFIX}_seeding_topology.txt 
python $MIGRATION ${PREFIX}_cp_output/${PREFIX}_all_results.txt > ${PREFIX}_cp_output/${PREFIX}_migration.txt

## ANALYSE SELECTION ON ORIGINAL AND MACHINA TOPOLOGY

python $SELECTION ${PREFIX}_cp_output/${PREFIX}_all_results.txt $ASV > ${PREFIX}_cp_output/${PREFIX}_selection.txt
python $SELECTION ${PREFIX}_cp_output/${PREFIX}_all_original_tree.txt $ASV > ${PREFIX}_cp_output/${PREFIX}_selection_original.txt
