#!/bin/bash

#source ~/miniconda3/etc/profile.d/conda.sh

run_evotracer=false
run_amplican=false
run_crispresso=false
run_machina=false
MIGRATION_MATRIX="NA"

if [[ $# -eq 0 ]] ; then
    ### Can uncomment below if the simulator.py is made to use two mutation rates
    #echo "Usage: sim_full_pipeline.sh --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> --migration <file_path>"
    echo "Usage: sim_full_pipeline.sh [--evotracer] [--amplican] [--crispresso2] [--machina] --out <out_name> --mutrate <float|comma-sep list of floats> --samples <int> [--migration <migration_matrix_filepath>]"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --evotracer) run_evotracer=true ;;
        --amplican) run_amplican=true ;;
        --crispresso2) run_crispresso=true ;;
        --machina) run_machina=true run_evotracer=true run_amplican=false run_crispresso=false;;
        -o|--out) NAME="$2"; shift ;;
        -m|--mutrate) MUTRATE="$2"; shift ;;
        -s|--samples) NUM_SAMPLES="$2"; shift ;;
        --migration) MIGRATION_MATRIX="$2"; shift ;;

    *) echo "Unknown parameter passed: $1"; echo "Usage: sim_full_pipeline.sh [--evo] [--amp] [--cresso] --out <out_name> --mutrate1 <float> --mutrate2 <float> --max-indel-size <int> --samples <int> [--migration <migration_matrix_filepath>]" ; exit 1 ;;
    esac
    shift
done

################################################################
### Run the simulator with the given command line parameters ###
################################################################

echo "Starting simulator..."
# you will likely want to create and activate a conda env using the provided yaml file
# conda env create -f env/simulate.yaml
conda activate simulate
eval "$(conda shell.bash hook)"

if [ ! -f "./scripts/simulator/sim_wrapper.sh" ]
then
    echo "Script ./scripts/simulator/sim_wrapper.sh not found. Exiting!"
    exit
fi

if [[ "$MIGRATION_MATRIX" != "NA" ]]; then
    ./scripts/simulator/sim_wrapper.sh --out ${NAME} --mutrate ${MUTRATE} --samples ${NUM_SAMPLES} --migration ${MIGRATION_MATRIX}
else
    ./scripts/simulator/sim_wrapper.sh --out ${NAME} --mutrate ${MUTRATE} --samples ${NUM_SAMPLES} --migration "NA"
fi

outputdir="sim_results_${NAME}/"
sim_dir="${outputdir}out_simulator_${NAME}/"

if [ -d "${sim_dir}" ]; then
  echo "Output directory ${sim_dir} already exists. Exiting!"
  exit
fi

mkdir ${sim_dir}

ls $outputdir| grep -v out_simulator_${NAME}| while read l; do mv $outputdir/${l} $outputdir/out_simulator_${NAME}/;done
#mv ${outputdir}/* ${sim_dir}
#echo "This message is okay"

if [[ "$MIGRATION_MATRIX" = "NA" ]]; then
    num_uniq_barcodes=$(python scripts/analysis/sim_barcode_count.py ${sim_dir}*tissue1*.fa)
else
    num_uniq_barcodes=$(python scripts/analysis/sim_barcode_count.py ${sim_dir}${NAME}.fa)
fi

sim_file="${sim_dir}${NAME}_mutations.tsv"

conda deactivate


# Run EvoTraceR package if flag is provided
if [[ "$run_evotracer" = true ]]; then
    echo "Starting EvoTraceR analysis..."
    #############################################
    ### Run EvoTraceR on the simulator output ###
    #############################################

    # you will likely want to create and activate a conda env using the provided yaml file
    # conda env create -f env/evotracer.yaml
    # then need to install R package in R from the library(devtools) for the EvoTraceR-parallelize repo and potentially fix "umi_tools" and "cassiopeia" issues if encountered
    conda activate evotracer

    if [ ! -f "./scripts/evotracer/run_evotracer_parallelize.R" ]
    then
        echo "Script ./scripts/evotracer/run_evotracer_parallelize.R not found. Exiting!"
        exit
    fi

    flash_path=$(which flash)
    if [ ! -x "$flash_path" ] ; then
        echo "FLASH executable not found in PATH. Exiting!"
        exit
    fi
    # ensure trimmomatic binary is executable
    trimmomatic_path=$(which trimmomatic-0.39.jar)
    if [ ! -x "$trimmomatic_path" ] ; then
        trimmomatic_path=$(which trimmomatic.jar)
        if [ ! -x "$trimmomatic_path" ] ; then
            echo "trimmomatic.jar executable not found in PATH. Exiting!"
            exit
        fi
    fi
    if [[ "$MIGRATION_MATRIX" = "NA" ]]; then
        evo_input_dir="./${sim_dir}all_samples_fastq"
    else
        evo_input_dir="./${sim_dir}tissue_specific_fastq"
    fi
    evo_output_dir="./${outputdir}out_evotracer_parallelize_${NAME}"

    #Rscript ./scripts/evotracer/run_evotracer_parallelize.R ${evo_input_dir} ${evo_output_dir} ${trimmomatic_path} ${flash_path} ${evotracer_path}
    Rscript ./scripts/evotracer/run_evotracer_parallelize.R ${evo_input_dir} ${evo_output_dir} ${trimmomatic_path} ${flash_path}

    num_uniq_mutations=$(wc -l ${sim_dir}${NAME}_mutations.tsv | awk '{print $1}')
    num_ASVs=$(tail -n +2 ${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv | cut -d ',' -f 1 | sort -u | wc -l)
    proportion=$(echo "scale=4;$num_ASVs/$num_uniq_barcodes" | bc)

    mr_write=${MUTRATE//,/ }

    conda deactivate

    evo_file="${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv"

    ############################################
    ### Calculate EvoTraceR alignment error of ASVs retrospectively to unique barcode match
    ############################################

    if [[ "$MIGRATION_MATRIX" = "NA" ]]; then
        conda activate evotracer
        ### Compare number of unique barcodes to number of collapsed ASVs
        echo "name,mutrate,num_samples,num_uniq_mutations,num_uniq_barcodes,num_ASVs,proportion_ASVs_barcodes,avg_alignment_error,num_barcodes_missing,proportion_barcodes_missing,proportion_barcodes_found,num_barcodes_found,proportion_ASVs_barcodes_found" >> ${outputdir}comparison_sim_evotracer_${NAME}.csv

        asv_filepath="${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv"
        sim_fasta_filepath="${sim_dir}MMUS1469_tissue1_BC10v0_MG_120419_${NAME}.fa"

        Rscript scripts/analysis/alignment_error_barcode_ASVs.R ${sim_fasta_filepath} ${asv_filepath} ./${outputdir}

        conda deactivate

        align_sum=$(awk -F',' 'NR>1{align_sum+=$2}END{print align_sum}' ${outputdir}alignment_errors.csv)
        avg_alignment_error=$(echo "scale=4;$align_sum/$num_ASVs" | bc)
        num_barcodes_missing=$(wc -l ${outputdir}barcode_no_ASV_alignment_sequences.txt | awk '{print $1}')
        proportion_barcodes_missing=$(echo "scale=4;$num_barcodes_missing/$num_uniq_barcodes" | bc)
        proportion_barcodes_found=$(echo "scale=4; 1 - $proportion_barcodes_missing" | bc)
        num_barcodes_found=$((num_uniq_barcodes - num_barcodes_missing))
        proportion_ASVs_barcodes_found=$(echo "scale=4;$num_ASVs/$num_barcodes_found" | bc)

        echo "${NAME},${mr_write},${NUM_SAMPLES},${num_uniq_mutations},${num_uniq_barcodes},${num_ASVs},${proportion},${avg_alignment_error},${num_barcodes_missing},${proportion_barcodes_missing},${proportion_barcodes_found},${num_barcodes_found},${proportion_ASVs_barcodes_found}" >> ${outputdir}comparison_sim_evotracer_${NAME}.csv
    fi
else 
    evo_file="NA"
fi

# Run MACHINA if flag is provided
if [[ "$run_machina" = true ]]; then
    #############################################
    ### Run Machina on the EvoTracer output ###
    #############################################

    echo "Starting Machina analysis..."
    # you will likely want to create and activate a conda env using the provided yaml file
    # conda env create -f env/machina.yaml
    conda activate machina

    if [ ! -f "./scripts/machina/run_machina.sh" ]
    then
        echo "Script ./scripts/machina/run_machina.sh not found. Exiting!"
        exit
    fi

    ASV="${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv"
    TREE="${evo_output_dir}/phylogeny_analysis/phylogeny_del_ins/tree_all_clones.newick"
    SPATH="./scripts/machina/"
    PREFIX="out_machina_${NAME}"
    RAW_TISSUES=$(head -n 1 ${MIGRATION_MATRIX})
    IFS=', ' read -r -a TISSUES <<< "${RAW_TISSUES}"
    PTISSUE="${TISSUES[1]}"

    timeout 3m ./scripts/machina/run_machina.sh --infile ${ASV} --tree ${TREE} --scripts ${SPATH} --prefix ${PREFIX} --primary-tissue ${PTISSUE} --keep-first-cp

    machina_dir="${outputdir}${PREFIX}/"
    mkdir ${machina_dir}
    rm -r ${PREFIX}_cp_output/data
    mv ${PREFIX}_cp_output/* ${machina_dir}/
    rm -r ${PREFIX}_cp_output/

    conda deactivate

    ##################################################
    ### Compare total true and inferred migrations ###
    ##################################################
    echo "Starting simulated ground truth and Machina comparison..."
    # Extract the specified column and count the number of True values for the true migrations total
    true_count=$(awk -F'\t' '{print $6}' ${sim_dir}${NAME}_tissues.tsv | grep -c True)

    # Extract numerical values from the "migrations" column and sum them for the inferred migrations total
    inferred_filepath="${machina_dir}${PREFIX}_migration.txt"
    if [ -e "$inferred_filepath" ]; then
        inferred_count=$(awk -F',' 'NR>1{sum+=$3} END{print sum}' ${inferred_filepath})
    else
        inferred_count=0
    fi
        
    #echo "Total inferred migrations in the ${NAME} simulation: $inferred_count"
    if (( true_count > 0 && inferred_count > 0 )); then
        proportion=$(echo "scale=4;$inferred_count/$true_count" | bc)
    elif (( true_count == 0 && inferred_count == 0 )); then
        proportion=1
    else
        proportion=0
    fi

    num_uniq_mutations=$(wc -l ${sim_dir}${NAME}_mutations.tsv | awk '{print $1}')

    avg_mut_age=$(awk '{ age += $NF } END { print age / NR }' ${sim_dir}${NAME}_mutations.tsv)

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
    done < <(tail -n +2 "${sim_dir}${NAME}_indel_character_matrix.tsv")
    avg_row_dropout=$(echo "scale=4; $total_dropout / $NUM_SAMPLES" | bc)
    avg_mut_sites_per_sample=$(echo "scale=4; $total_mut_sites / $NUM_SAMPLES" | bc)
    avg_proportion_mut_sites=$(echo "scale=4; $avg_mut_sites_per_sample / $num_sites" | bc)

    echo "${NAME},${mr_write},${NUM_SAMPLES},${MIGRATION_MATRIX},${num_uniq_mutations},${uniq_mut_per_site},${total_mut_sites},${avg_mut_sites_per_sample},${avg_proportion_mut_sites},${total_dropout},${avg_row_dropout},${avg_mut_age},${true_count},${inferred_count},${proportion}" >> ${outputdir}comparison_inferred_true_migration_${NAME}.csv

    conda activate simulate

    python ./scripts/analysis/compare_migrations_simtrue_machina.py "${sim_dir}${NAME}_tissues.tsv" ${inferred_filepath} "${NAME}" "${outputdir}" ${MIGRATION_MATRIX}

    conda deactivate

fi


# Run ampliCan package if flag is provided
if [[ "$run_amplican" = true ]]; then
    ### ampliCan can only run if the passed mutation rate vector has all 0 values except for in position 5 out of 10 such as [0,0,0,0,0.1,0,0,0,0,0]
    ### other positions can be setup for analysis with the limitation of necessitating perfect guide match (not experimental imperfect guide) and false primer design

    echo "Starting ampliCan analysis..."
    #############################################
    ### Run ampliCan on the simulator output ###
    #############################################
    conda activate amplican_env

    amplican_dir="${outputdir}out_amplican_${NAME}/"
    mkdir ${amplican_dir}

    echo "ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor" >> ${amplican_dir}amplican_config_${NAME}.csv

    ## fake primers to use amplican on simulated fastq files for binding sites of primers and guideRNA to site 5
    guideRNA="TTTACATACTCGTTCAACCG"
    forward_primer="TCTACACGCGCGTTCAAC"
    reverse_primer="GATCCACGGGTGAATGCA"
    amplicon="tctacacgcgcgttcaaccgaggaaaactacacacacgttcaaccacggttttttacacacgcattcaaccacggactgctacacacgcactcaaccgtggataTTTACATACTCGTTCAACCGtggattgttacacccgcgttcaaccagggtcagatacacccacgttcaaccgtggtactatactcgggcattcaaccgcggctttctgcacacgcctacaaccgcggaactatacacgtgcattcacccgtggatc"

    echo "tissue1,barcode_1,MMUS1469_tissue1_BC10v0_MG_120419_${NAME}_R1.fastq,MMUS1469_tissue1_BC10v0_MG_120419_${NAME}_R2.fastq,group1,1,${guideRNA},${forward_primer},${reverse_primer},0,${amplicon}," >> ${amplican_dir}amplican_config_${NAME}.csv

    timeout 20m Rscript ./scripts/amplican/run_amplican.R ${sim_dir}all_samples_fastq/ ${amplican_dir}amplican_config_${NAME}.csv ${amplican_dir}
    
    conda deactivate

    amp_file="${amplican_dir}alignments/events_filtered_shifted_normalized.csv"
    amp_count_file="${amplican_dir}barcode_reads_filters.csv"
else 
    amp_file="NA"
    amp_count_file="NA"
fi


# Run CRISPResso package if flag is provided
if [[ "$run_crispresso" = true ]]; then

    echo "Starting CRISPResso2 analysis..."
    #############################################
    ### Run CRISPResso2 on the simulator output ###
    #############################################

    conda activate crispresso2_env

    cresso_dir="${outputdir}out_CRISPResso2_${NAME}/"
    mkdir ${cresso_dir}

    fastq_r1="${sim_dir}all_samples_fastq/MMUS1469_tissue1_BC10v0_MG_120419_${NAME}_R1.fastq"
    fastq_r2="${sim_dir}all_samples_fastq/MMUS1469_tissue1_BC10v0_MG_120419_${NAME}_R2.fastq"
    amplicon="TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC"
    guideRNA="TTTACATACTCGTTCAACCG" ## fake guide for site 5 comparison

    CRISPResso --ignore_substitutions --suppress_report --suppress_plots --fastq_r1 ${fastq_r1} --fastq_r1 ${fastq_r2} --amplicon_seq ${amplicon} --guide_seq ${guideRNA} -o ${cresso_dir}

    cresso_file="${cresso_dir}CRISPResso_on_MMUS1469_tissue1_BC10v0_MG_120419_${NAME}_R2/Alleles_frequency_table_around_sgRNA_TTTACATACTCGTTCAACCG.txt"
    conda deactivate
else
    cresso_file="NA"
fi

if [[ "$run_evotracer" = true ]] && [[ "$run_amplican" = true ]] && [[ "$run_crispresso" = true ]]; then
    conda activate simulate
    python ./scripts/analysis/compare_amplican_evotracer_crispresso2.py ${NAME} "${mr_write}" ${NUM_SAMPLES} ${sim_file} ${amp_file} ${outputdir} ${amp_count_file} ${evo_file} ${cresso_file}
                
    conda deactivate
fi

 
