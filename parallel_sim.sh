#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh


parallel_sim_name="explore_parameters"
# Set the migration matrix values that will be tested
rare_mm="data/rare_migration_prob_matrix.csv"
equal_mm="data/equal_migration_prob_matrix.csv"
true_mm="data/true_migration_prob_matrix.csv"

mm_array=(${rare_mm} ${equal_mm} ${true_mm})

# Set the mutation rates to explore
mr_array=(0.01 0.025 0.05 0.1 0.2)

# Set the sample sizes to explore
ss_array=(100 1000 10000 50000 100000)

# Create an array to store all the commands
commands=()

# Loop over all combinations of the three arrays and add the command to the commands array
for mm in ${mm_array[@]}
do
  for mr in ${mr_array[@]}
  do
    for ss in ${ss_array[@]}
    do
      # Construct the command to submit with ParaFly
      simmid="sim_${mr}_${ss}_${mm}"
      cmd="./sim_full_pipeline.sh --out ${simmid} --mutrate ${mr} --max-indel-size 5 --samples ${ss} --migration ${mm}"

      # Add the command to the commands array
      commands+=("$cmd")
    done
  done
done
echo "There are ${#commands[@]} commands to be submitted."

# Submit the commands in batches of 20 using ParaFly
batch_size=20
num_batches=$(( (${#commands[@]} + batch_size - 1) / batch_size ))

for ((i=0; i<$num_batches; i++))
do
  # Calculate the range of commands to submit in this batch
  start_index=$((i * batch_size))
  end_index=$((start_index + batch_size - 1))

  # Slice the commands array to get the commands for this batch
  batch_commands=("${commands[@]:$start_index:$batch_size}")

  # Construct the ParaFly command to submit this batch of commands with 20 CPUs
  cmd="ParaFly -c '${batch_commands[*]}' -CPU 20"

  # Print the ParaFly command to the console
  echo "Submitting batch $((i+1)) of $num_batches: $cmd"

  # Submit the batch with ParaFly
  eval $cmd
done

# Create the main output file with the header
echo "name,mutrate,max_indel_size,num_samples,migration_matrix,true_migrations,inferred_migrations,proportion" > ${parallel_sim_name}_all_output.csv

# Append the second line of each csv file to the main file
for f in sim_results*/comparison*; do tail -n +2 $f >> ${parallel_sim_name}_all_output.csv; done

# Remove all the individual csv files
rm -r sim_results*

