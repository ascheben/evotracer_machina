ASV="$1"
TREE="$2"

TRAV="/local/storage/no-backup/scheben-scratch/proca/machina/traverse_split.py"
GET="/local/storage/no-backup/scheben-scratch/proca/machina/get_results.sh"
MACHINA="~/miniconda3/envs/machina/bin/pmh_tr"


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

# To Do: 
# 1) all or most of the analysis scripts use a loop that ignores the final CP in the input results
# 2) produce CSV formatted output from analysis scripts

python /local/storage/no-backup/scheben-scratch/proca/machina/print_seeding_topology.py all_results_migration_optimal_split_22082022.txt | grep seeding| sed 's/ .*{/ /'| sed 's/}//'| cut -d' ' -f1,4,7,10,12,14,17 | sed 's/,//g' | tr ' ' '\t'

python /local/storage/no-backup/scheben-scratch/proca/machina/selection_test.py all_results.txt asv_stat.csv

python /local/storage/no-backup/scheben-scratch/proca/machina/count_migrations.py all_results.txt

