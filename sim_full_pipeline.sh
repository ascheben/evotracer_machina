#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

################################################################
### Run the simulator with the given command line parameters ###
################################################################

echo "Starting simulator..."
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
conda activate simulate

if [[ $# -eq 0 ]] ; then
    ### Can uncomment below if the simulator.py is made to use two mutation rates
    #echo "Usage: sim_full_pipeline.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> --migration <file_path>"
    echo "Usage: sim_full_pipeline.sh --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> --migration <file_path>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -o|--out) NAME="$2"; shift ;;
        -m|--mutrate) MUTRATE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;
        -mm|--migration) MIGRATION_MATRIX="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_full_pipeline.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> --migration <file_path>" ; exit 1 ;;
    esac
    shift
done

if [ ! -f "./sim_wrapper.sh" ]
then
    echo "Script ./sim_wrapper.sh not found. Exiting!"
    exit
fi

./sim_wrapper.sh --out ${NAME} --mutrate ${MUTRATE} --samples ${NUM_SAMPLES} --migration ${MIGRATION_MATRIX}

outputdir="sim_results_${NAME}/"
sim_dir="${outputdir}out_simulator_${NAME}/"
mkdir ${sim_dir}
mv ${outputdir}* ${sim_dir}
echo "This message is okay"

conda deactivate

#############################################
### Run Evotracer on the simulator output ###
#############################################

echo "Starting EvoTraceR analysis..."
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/r_env.yaml
# then need to install R package in R from the library(devtools) for the EvoTraceR-parallelize repo and potentially fix "umi_tools" and "cassiopeia" issues if encountered

conda activate r_env

if [ ! -f "./scripts/run_evotracer_parallelize.R" ]
then
    echo "Script ./scripts/run_evotracer_parallelize.R not found. Exiting!"
    exit
fi

# Need to manually specify paths here for the following variables
trimmomatic_path="/home/staklins/miniconda3/envs/r_env/share/trimmomatic-0.39-2/trimmomatic.jar"
flash_path="/home/staklins/miniconda3/envs/r_env/bin/flash"
evotracer_path="../EvoTraceR-parallelize"

evo_input_dir="./${sim_dir}tissue_specific_fastq"
evo_output_dir="./${outputdir}out_evotracer_parallelize_${NAME}"

Rscript ./scripts/run_evotracer_parallelize.R ${evo_input_dir} ${evo_output_dir} ${trimmomatic_path} ${flash_path} ${evotracer_path}

conda deactivate

#############################################
### Run Machina on the EvoTracer output ###
#############################################

echo "Starting Machina analysis..."
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/machina.yaml
conda activate machina

if [ ! -f "./run_pipeline.sh" ]
then
    echo "Script ./run_pipeline.sh not found. Exiting!"
    exit
fi

ASV="${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv"
TREE="${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick"
SPATH="./scripts/"
PREFIX="out_machina_${NAME}"
RAW_TISSUES=$(head -n 1 ${MIGRATION_MATRIX})
IFS=', ' read -r -a TISSUES <<< "${RAW_TISSUES}"
PTISSUE="${TISSUES[1]}"

timeout 3m ./run_pipeline.sh --infile ${ASV} --tree ${TREE} --scripts ${SPATH} --prefix ${PREFIX} --primary-tissue ${PTISSUE} --keep-first-cp

mkdir ${outputdir}${PREFIX}
mv ${PREFIX}_cp_output/* ${outputdir}${PREFIX}/
rm -r ${PREFIX}_cp_output/

conda deactivate


### Retain only migration related files and delete the rest
# Prune files for simulator output
mv ${sim_dir}${NAME}_true_tissues.nwk ${outputdir}
mv ${sim_dir}${NAME}_tissues.tsv ${outputdir}
mv ${sim_dir}${NAME}_mutations.tsv ${outputdir}
mv ${sim_dir}${NAME}_indel_character_matrix.tsv ${outputdir}
rm -r ${sim_dir}

# Prune all EvoTraceR output files
rm -r ${evo_output_dir}

# Prune Machina intermediate output files
mv ${outputdir}${PREFIX}/${PREFIX}* ${outputdir}
rm -r ${outputdir}${PREFIX}

##################################################
### Compare total true and inferred migrations ###
##################################################
echo "Starting simulated ground truth and Machina comparison..."
# Extract the specified column and count the number of True values for the true migrations total
true_count=$(awk -F'\t' '{print $5}' ${outputdir}${NAME}_tissues.tsv | grep -c True)
#echo "Total true migrations in the ${NAME} simulation: $true_count"

# Extract numerical values from the "migrations" column and sum them for the inferred migrations total
inferred_count=$(awk -F',' 'NR>1{sum+=$3} END{print sum}' ${outputdir}${PREFIX}_migration.txt)
#echo "Total inferred migrations in the ${NAME} simulation: $inferred_count"
proportion=$(echo "scale=4;$inferred_count/$true_count" | bc)

num_uniq_mutations=$(wc -l ${outputdir}${NAME}_mutations.tsv | awk '{print $1}')

avg_mut_age=$(awk '{ age += $NF } END { print age / NR }' ${outputdir}${NAME}_mutations.tsv)

# Write both true and inferred to an output csv with the input parameters stored
echo "name,mutrate,num_samples,migration_matrix,num_uniq_mutations,uniq_mut_per_site,total_mut_sites,avg_mut_sites_per_sample,avg_proportion_mut_sites,total_dropout,dropout_per_sample,avg_mutation_age,true_migrations,inferred_migrations,proportion" >> ${outputdir}comparison_inferred_true_migration_${NAME}.csv
mr_write=${MUTRATE//,/ }
# Loop over the array and count the non-zero values
IFS=' ' read -r -a mut_array <<< "$mr_write"
num_sites=0
for i in "${mut_array[@]}"
do
  if [[ $i != '0' ]]
  then
    num_sites=$((num_sites+1))
  fi   
done
uniq_mut_per_site=$(echo "scale=4; $num_uniq_mutations / $num_sites" | bc)

# Initialize a variable to keep track of the total number of -1 values
total_dropout=0
declare -a dropout_array
total_mut_sites=0
declare -a mut_array

while read line; do
    # Count the number of -1 values in the row
    row_dropout=$(echo "$line" | cut -f 2- | tr '\t' '\n' | grep -e "^-1$" | wc -l)
    total_dropout=$((total_dropout + row_dropout))
    dropout_array+=("$row_dropout")
    # Count non 0 and non -1 to get number of mutated sites per row
    row_mut_count=$(echo "$line" | cut -f 2- | tr '\t' '\n' | grep -v -e "^$" -e "^-1$" -e "^0$" | wc -l)
    total_mut_sites=$((total_mut_sites + row_mut_count))
    mut_array+=("$row_mut_count")
done < <(tail -n +2 "${outputdir}${NAME}_indel_character_matrix.tsv")
avg_row_dropout=$(echo "scale=4; $total_dropout / $NUM_SAMPLES" | bc)
avg_mut_sites_per_sample=$(echo "scale=4; $total_mut_sites / $NUM_SAMPLES" | bc)
avg_proportion_mut_sites=$(echo "scale=4; $avg_mut_sites_per_sample / $num_sites" | bc)

echo "${NAME},${mr_write},${NUM_SAMPLES},${MIGRATION_MATRIX},${num_uniq_mutations},${uniq_mut_per_site},${total_mut_sites},${avg_mut_sites_per_sample},${avg_proportion_mut_sites},${total_dropout},${avg_row_dropout},${avg_mut_age},${true_count},${inferred_count},${proportion}" >> ${outputdir}comparison_inferred_true_migration_${NAME}.csv

conda activate simulate

python ./scripts/compare_migrations_simtrue_machina.py "${outputdir}${NAME}_tissues.tsv" "${outputdir}${PREFIX}_migration.txt" "${NAME}" "${outputdir}"

conda deactivate

rm ${outputdir}${NAME}*
rm ${outputdir}${PREFIX}*

