#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

parallel_sim_name="simEvoCompareBulkMutrates_FlexibleTP1bpEndMismatch_7_17_23"
running_machina=false   ### SET TO TRUE IF RUNNING MACHINA

# Set the mutation rates to explore

# Empirical mutation rates
#mr1=(0.3173305,0.0009823,0.00000008,0.00000162,0.00000335,0.18609872,0.00000254,0,0,0)

# Easy alignment mode mutation rates
#mr1=(0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264)
#mr1=(0,0,0,0,0.1,0,0,0,0,0)
#mr_array=("$mr1")
# mr1=(0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001)
# mr2=(0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005)
mr1=(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
mr2=(0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02)
mr3=(0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03)
mr4=(0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04)
mr5=(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
mr6=(0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06)
mr7=(0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07)
mr8=(0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08)
mr9=(0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09)
mr10=(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
mr11=(0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15)
mr12=(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
mr13=(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3)
mr14=(0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4)
mr15=(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
# mr18=(0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
# mr19=(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7)
# mr20=(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)
# mr21=(0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9)
# mr22=(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
mr_array=("$mr1" "$mr2" "$mr3" "$mr4" "$mr5" "$mr6" "$mr7" "$mr8" "$mr9" "$mr10" "$mr11" "$mr12" "$mr13" "$mr14" "$mr15")

# Set the sample sizes to explore
ss_array=(100)

# Set the migration matrix values that will be tested
# ultra_rare_mm="data/ultrarare_migration_prob_matrix.csv"
# rare_mm="data/rare_migration_prob_matrix.csv"
#moderate_mm="data/moderate_migration_prob_matrix.csv"
# true_mm="data/true_migration_prob_matrix.csv"
# high_mm="data/high_migration_prob_matrix.csv"
#equal_mm="data/equal_migration_prob_matrix.csv"
# mm_array=(${rare_mm} ${true_mm} ${high_mm})
# mm_array=(${true_mm})

mm_1="data/migration_matrices/migration_matrix_001.csv"  ## 0.001 migration rate split between all tissues
mm_2="data/migration_matrices/migration_matrix_005.csv"   ## 0.005 migration rate split between all tissues
mm_3="data/migration_matrices/migration_matrix_01.csv"   ## 0.01 migration rate split between all tissues
mm_4="data/migration_matrices/migration_matrix_025.csv"   ## 0.025 migration rate split between all tissues
mm_5="data/migration_matrices/migration_matrix_05.csv"   ## 0.05 migration rate split between all tissues
mm_6="data/migration_matrices/migration_matrix_075.csv"   ## 0.075 migration rate split between all tissues
mm_7="data/migration_matrices/migration_matrix_1.csv"   ## 0.1 migration rate split between all tissues
mm_8="data/migration_matrices/migration_matrix_25.csv"   ## 0.25 migration rate split between all tissues
mm_9="data/migration_matrices/migration_matrix_50.csv"   ## 0.50 migration rate split between all tissues
mm_10="data/migration_matrices/migration_matrix_75.csv"   ## 0.75 migration rate split between all tissues
mm_11="data/migration_matrices/migration_matrix_99.csv"   ## 0.99 migration rate split between all tissues
mm_array=(${mm_1} ${mm_2} ${mm_3} ${mm_4} ${mm_5} ${mm_6} ${mm_7} ${mm_8} ${mm_9} ${mm_10} ${mm_11})



if [[ "$running_machina" = false ]]; then
  # Setup headers for recording the input parameters in a csv
  mr_header="mutation_rate"
  ss_header="sample_size"

  # Determine the maximum length of the arrays for input parameter csv writing
  max_len=${#mm_array[@]}
  if [ ${#mr_array[@]} -gt $max_len ]; then
    max_len=${#mr_array[@]}
  fi
  if [ ${#ss_array[@]} -gt $max_len ]; then
    max_len=${#ss_array[@]}
  fi

  # Write the array values to a CSV file
  echo "${mr_header},${ss_header}" > input_ranges_${parallel_sim_name}.csv

  for ((i=0; i<$max_len; i++)); 
    do
      mr=${mr_array[$i]:-}
      mr=${mr//,/ }
      ss=${ss_array[$i]:-}
      echo "${mr},${ss}" >> input_ranges_${parallel_sim_name}.csv
    done

  # Create an array to store all the commands
  commands=()

  echo "Setting up all parameter combinations..."
  # Loop over all combinations of the three arrays and add the command to the commands array
  j=0

  for mr in ${mr_array[@]}
  do
    for ss in ${ss_array[@]}
    do
      for ((rep=0; rep<100; rep++))   ### Can set the amount of repeat simulations with the same parameters
      do
      simmid="sim${j}"
      cmd="./sim_full_pipeline.sh --evotracer --out ${simmid} --mutrate ${mr} --samples ${ss}"
      commands+=("$cmd")
      j=$((j+1))
      done
    done
  done
fi



if [[ "$running_machina" = true ]]; then
  # Setup headers for recording the input parameters in a csv
  mr_header="mutation_rate"
  ss_header="sample_size"
  mm_header="migration_matrix"

  # Determine the maximum length of the arrays for input parameter csv writing
  max_len=${#mm_array[@]}
  if [ ${#mr_array[@]} -gt $max_len ]; then
    max_len=${#mr_array[@]}
  fi
  if [ ${#ss_array[@]} -gt $max_len ]; then
    max_len=${#ss_array[@]}
  fi

  # Write the array values to a CSV file
  echo "${mr_header},${ss_header},${mm_header}" > input_ranges_${parallel_sim_name}.csv

  for ((i=0; i<$max_len; i++)); 
    do
      mr=${mr_array[$i]:-}
      mr=${mr//,/ }
      ss=${ss_array[$i]:-}
      mm=${mm_array[$i]:-}
      echo "${mr},${ss},${mm}" >> input_ranges_${parallel_sim_name}.csv
    done

  # Create an array to store all the commands
  commands=()

  echo "Setting up all parameter combinations..."
  # Loop over all combinations of the three arrays and add the command to the commands array
  j=0
  for mm in ${mm_array[@]}
  do
    for mr in ${mr_array[@]}
    do
      for ss in ${ss_array[@]}
      do
        for ((rep=0; rep<100; rep++))   ### Can set the amount of repeat simulations with the same parameters
        do
        simmid="sim${j}"
        cmd="./sim_full_pipeline.sh --machina --out ${simmid} --mutrate ${mr} --samples ${ss} --migration ${mm}"
        commands+=("$cmd")
        j=$((j+1))
        done
      done
    done
  done
fi

echo "There are ${#commands[@]} commands to be submitted."

# Create the main output files with the headers. These can be commented out if not running specific packages in the sim_full_pipeline.sh command line optional flag inputs.
echo "name,mutrate,num_samples,num_uniq_mutations,num_uniq_barcodes,num_ASVs,proportion_ASVs_barcodes,avg_alignment_error,num_barcodes_missing,proportion_barcodes_missing,proportion_barcodes_found,num_barcodes_found,proportion_ASVs_barcodes_found" > output_sim_evotracer_${parallel_sim_name}.csv
echo "sim_name,mutrate,num_samples,sim_mutations,count_sim_mutations,evotracer_mutations,count_evotracer_mutations,evo_precision,evo_recall,evo_f1,amplican_mutations,count_amplican_mutations,amp_precision,amp_recall,amp_f1,crispresso2_mutations,count_crispresso2_mutations,crispresso2_precision,crispresso2_recall,crispresso2_f1" > output_simEvoAmpCressoCompare_${parallel_sim_name}.csv
if [["$running_machina" = true]]; then
  echo "name,mutrate,num_samples,migration_matrix,num_uniq_mutations,uniq_mut_per_site,total_mut_sites,avg_mut_sites_per_sample,avg_proportion_mut_sites,total_dropout,dropout_per_sample,avg_mutation_age,true_migrations,inferred_migrations,proportion" > output_machina_all_${parallel_sim_name}.csv
  echo "machina_sim_name,migration_matrix,true_positives,false_positives,false_negatives,precision,recall,f1_score" > output_statistics_machina_all_${parallel_sim_name}.csv
fi

par_results="${parallel_sim_name}_parallel_sim_results"
mkdir ${par_results}
par_results_data="${par_results}/data"
mkdir ${par_results_data}

# Submit the commands in batches of 20 using ParaFly
batch_size=20
num_batches=$(((${#commands[@]} / $batch_size) + 1))

echo "Beginning submission of all combinations in parallel batches of ${batch_size}..."
conda activate machina
for ((i=0; i<$num_batches; i++))
do
  echo "Submitting batch $((i+1))/$((num_batches)):"
  # Calculate the range of commands to submit in this batch
  start_index=$((i * batch_size))
  end_index=$((start_index + batch_size - 1))

  # Check if the end index is greater than or equal to the length of the commands array
  if [ $end_index -ge ${#commands[@]} ]; then
    batch_commands=("${commands[@]:${start_index}}") # Slice from start index to end of commands array
  else
    batch_commands=("${commands[@]:${start_index}:${batch_size}}") # Slice from start index to batch_size
  fi

  for command in "${batch_commands[@]}"
  do
    echo "${command}" >> "${i}.cmd"
    echo "${command}"
  done
  # Print the ParaFly command to the console
  ParaFly -CPU ${batch_size} -c ${i}.cmd
  rm ${i}.cmd
  rm ${i}.cmd.completed
  # Append the second line of each csv file to the main files
  for f in sim_results*/comparison_sim*; do tail -n +2 $f >> output_sim_evotracer_${parallel_sim_name}.csv; done
  for f in sim_results*/compare_*; do tail -n +2 $f >> output_simEvoAmpCressoCompare_${parallel_sim_name}.csv; done

  if [["$running_machina" = true]]; then
    for f in sim_results*/comparison_inferred*; do tail -n +2 $f >> output_machina_all_${parallel_sim_name}.csv; done
    for f in sim_results*/statistics*; do tail -n +2 $f >> output_statistics_machina_all_${parallel_sim_name}.csv; done
  fi

  mv sim_results_sim* ${par_results_data}/
  #rm -r sim_results*
done

echo "All batches are done."

conda deactivate

mv *${parallel_sim_name}.csv ${par_results}/


### Use below to print all proportions data and remove sim data for use with scipy optimize
### Remove this if not using scipy since this deletes all output and prints to the terminal only what scipy needs
#proportions_scipy=$(tail -n +2 ${par_results}/output_all_${parallel_sim_name}.csv | cut -d ',' -f 15 | sed ':a;N;$!ba;s/\n/,/g')
#echo "proportions = [$proportions_scipy]"
#rm -r ${par_results}/


