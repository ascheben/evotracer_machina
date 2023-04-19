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
    echo "Usage: sim_full_pipeline.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> --migration <file_path>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -o|--out) NAME="$2"; shift ;;
        -m1|--mutrate1) MUTRATE1="$2"; shift ;;
        -m2|--mutrate2) MUTRATE2="$2"; shift ;;
        -i|--max-indel-size) MAX_INDEL_SIZE="$2"; shift ;;
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

./sim_wrapper.sh --out ${NAME} --mutrate1 ${MUTRATE1} --mutrate2 ${MUTRATE2} --max-indel-size ${MAX_INDEL_SIZE} --samples ${NUM_SAMPLES} --migration ${MIGRATION_MATRIX}

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

./run_pipeline.sh --infile ${ASV} --tree ${TREE} --scripts ${SPATH} --prefix ${PREFIX} --primary-tissue ${PTISSUE}

mkdir ${outputdir}${PREFIX}
mv ${PREFIX}_cp_output/* ${outputdir}${PREFIX}/
rm -r ${PREFIX}_cp_output/

conda deactivate


### Retain only migration related files and delete the rest
# Prune files for simulator output
mv ${sim_dir} ${NAME}_true_tissues.nwk ../
mv ${sim_dir} ${NAME}_tissues.tsv ../
rm -r ${sim_dir}

# Prune all EvoTraceR output files
rm -r ${evo_output_dir}

# Prune Machina intermediate output files
mv ${outputdir}${PREFIX}/${PREFIX}* ..
rm -r ${outputdir}${PREFIX}


