#!/bin/bash

absolute_path () {
 abspath=`echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"`
 echo $abspath
}


if [[ $# -eq 0 ]] ; then
    echo "Usage: run_pipeline.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts> --prefix <output_prefix> --primary-tissue <tissue> [--keep-first-cp]"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--infile) ASV="$2"; shift ;;
        -t|--tree) TREE="$2"; shift ;;
        -s|--scripts) SPATH="$2"; shift ;;
        -p|--prefix) PREFIX="$2"; shift ;;
        -o|--primary-tissue) PTISSUE="$2"; shift ;;
        -c|--cutoff) CUTOFF="$2"; shift ;;
        -k|--keep-first-cp) KEEP=true;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: run_pipeline.sh --infile <asv_stats> --tree <newick_tree> --scripts </path/to/scripts> --primary-tissue <tissue> [--keep-first-cp]" ; exit 1 ;;
    esac
    shift
done

if ! command -v pmh_tr &> /dev/null
then
    echo "pmh_tr from the MACHINA installation could not be found. Exiting!"
    exit
fi

if [ -z "${CUTOFF}" ]
then
    CUTOFF=1
fi

if [ -z "${ASV}" ]
then
    echo "Missing --infile parameter. Exiting!"
    exit
fi
if [ -z "${TREE}" ]
then
     echo "Missing --tree parameter. Exiting!"
    exit
fi
if [ -z "${SPATH}" ]
then
    echo "Missing --scripts parameter. Exiting!"
    exit
fi
if [ -z "${PREFIX}" ]
then
    echo "Missing --prefix parameter. Exiting!"
    exit
fi
if [ -z "${PTISSUE}" ]
then
    echo "Missing --primary-tissue parameter. Exiting!"
    exit
fi

# PATHS TO SCRIPTS 
THREADS=1
BATCHES=1
BIG_CP_THRESHOLD="${CUTOFF//[$'\t\r\n ']}"

SPATH=`absolute_path "${SPATH}"` 
TREE=`absolute_path "${TREE}"`
ASV=`absolute_path "${ASV}"`
TRAV="${SPATH}/traverse_split.py"
GET="${SPATH}/get_results.sh"
GETOLD="${SPATH}/get_results_old2new.sh"
MACHINA="pmh_tr"
TOPOLOGY="${SPATH}/print_seeding_topology.py"
SELECTION="${SPATH}/selection_tree_test.py"
MIGRATION="${SPATH}/count_migrations.py"
ADD_INFO="${SPATH}/add_freq_prob_to_results.py"

## PREPROCESS INPUT DATA ##

if [[ -n "${BIG_CP_THRESHOLD//[0-9]/}" ]]; then
    echo "Value for cutoff parameter is not an integer!"
    exit
fi

if [ -d "${PREFIX}_cp_output" ]; then
  echo "Output directory ${PREFIX}_cp_output already exist. Exiting!"
  exit
fi

mkdir ${PREFIX}_cp_output
cd ${PREFIX}_cp_output
# Extract key ASV columns
#asv_names,sample,group
cut -d',' -f1,2,30 ${ASV} | sed '/^$/d' > ${PREFIX}_asv_sample_group.csv
# exclude miscelleneaous group CP00
if [ "$KEEP" = true ] ; then
    remove="^$"
else
    first_cp =`cut -f3 -d',' ${PREFIX}_asv_sample_group.csv|tail -n +2| head -1`
    remove="^${first_cp}$"
fi
cut -f3 -d',' ${PREFIX}_asv_sample_group.csv|tail -n +2| sort| uniq| egrep -v "$remove"| sed '/^$/d' > ${PREFIX}_CP_list.txt
# prepare machina input files for each CP
touch ${PREFIX}_big_CP_list.txt
while read l; do
    ${TRAV} ${TREE} ${PREFIX}_asv_sample_group.csv ${l} ${PTISSUE}
    num_labels=`sed -n '$=' ${l}_labels_split.txt`
    if [[ $num_labels -gt $BIG_CP_THRESHOLD ]]; then
    #if [ "$num_labels" -gt 30 ]; then
        echo ${l} >> ${PREFIX}_big_CP_list.txt
    fi
done<${PREFIX}_CP_list.txt
# remove CPs that failed traversal due to unexpected clade topology
#grep -v -f FailedCP.txt ${PREFIX}_CP_list.txt > ${PREFIX}_CP_list.tmp
#mv ${PREFIX}_CP_list.tmp ${PREFIX}_CP_list.txt
# make output directory for each CP
while read l; do mkdir ${l}_split;done<${PREFIX}_CP_list.txt
   
# prepare input for selection scan on original tree
#mkdir ${PREFIX}_cp_output
for l in *_tree_split.txt; do cat $l |sed "s/^/${l} tree /"| sed 's/_tree_split.txt//';done | grep -v CP00 > ${PREFIX}_all_original_tree.txt

## RUN MACHINA ##     
# make command file for machina
# speed up MACHINA a bit by adding "-m 3"
grep -f ${PREFIX}_big_CP_list.txt ${PREFIX}_CP_list.txt| while read l; do echo "${MACHINA} -OLD -t ${THREADS} -m 3 -o ${l}_split -c ${l}_colors.txt -p ${PTISSUE} ${l}_tree_split.txt ${l}_labels_split.txt &> ${l}_split/results.txt";done >> ${PREFIX}_machina.cmd
grep -v -f ${PREFIX}_big_CP_list.txt ${PREFIX}_CP_list.txt| while read l; do echo "${MACHINA} -t ${THREADS} -m 3 -o ${l}_split -c ${l}_colors.txt -p ${PTISSUE} ${l}_tree_split.txt ${l}_labels_split.txt &> ${l}_split/results.txt";done >> ${PREFIX}_machina.cmd

# run machina in parallel
ParaFly -CPU ${BATCHES} -c ${PREFIX}_machina.cmd

# parse results from each machina output dir
grep -f ${PREFIX}_big_CP_list.txt ${PREFIX}_CP_list.txt| while read l; do ${GETOLD} $l ${PTISSUE} ${SPATH};done | tr '\t' ' '>> ${PREFIX}_all_results.txt
grep -v -f ${PREFIX}_big_CP_list.txt ${PREFIX}_CP_list.txt| while read l; do ${GET} $l ${PTISSUE} ${SPATH};done | tr '\t' ' '>> ${PREFIX}_all_results.txt

# move intermediate output
#mv CP* ${PREFIX}_cp_output
#mv ${PREFIX}_big_CP_list.txt ${PREFIX}_cp_output
#mv ${PREFIX}_machina* ${PREFIX}_cp_output
#mv ${PREFIX}_asv_sample_group.csv ${PREFIX}_cp_output
#mv ${PREFIX}_CP_list.txt ${PREFIX}_cp_output
## ANALYSE INFERRED TOPOLOGY

python $TOPOLOGY ${PREFIX}_all_results.txt ${PTISSUE} > ${PREFIX}_seeding_topology.txt 
python $MIGRATION ${PREFIX}_all_results.txt > ${PREFIX}_migration.txt

## ANALYSE SELECTION ON ORIGINAL AND MACHINA TOPOLOGY

#python $SELECTION ${PREFIX}_cp_output/${PREFIX}_all_results.txt $ASV > ${PREFIX}_cp_output/${PREFIX}_selection.txt
python $SELECTION ${PREFIX}_all_original_tree.txt $ASV| grep "^test" > ${PREFIX}_selection_original_test.txt
python $SELECTION ${PREFIX}_all_original_tree.txt $ASV| grep "^expansion" > ${PREFIX}_selection_original_expansion.txt
python $ADD_INFO ${PREFIX}_migration.txt $ASV ${PREFIX}_all_results.txt > ${PREFIX}_all_results_extended.txt

# Clean up
mkdir data
mv CP* data/
mv *list.txt data/
mv *cmd* data/
cd ..
