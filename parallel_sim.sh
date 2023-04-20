#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

parallel_sim_name="explore_parameters"
# Set the migration matrix values that will be tested
rare_mm="data/rare_migration_prob_matrix.csv"
equal_mm="data/equal_migration_prob_matrix.csv"
true_mm="data/true_migration_prob_matrix.csv"

# Set the mutation matrix filepaths
#mm_array=(${rare_mm} ${equal_mm} ${true_mm})
mm_array=(${equal_mm} ${true_mm})

# Set the mutation rates to explore
mr_array=(0.01 0.025 0.05 0.1 0.2)

# Set the sample sizes to explore
#ss_array=(100 1000 10000 50000 100000)
ss_array=(100 200)

# Create an array to store all the commands
commands=()

echo "Setting up all parameter combinations..."
# Loop over all combinations of the three arrays and add the command to the commands array
for mm in ${mm_array[@]}
do
  for mr in ${mr_array[@]}
  do
    for ss in ${ss_array[@]}
    do
      # Construct the command to submit with ParaFly
      mm_type=$(echo $mm | awk -F/ '{print $2}' | awk -F_ '{print $1}')
      mr_value=$(echo $mr | awk -F. '{print $2}')
      simmid="sim_${mr_value}_${ss}_${mm_type}"
      cmd="./sim_full_pipeline.sh --out ${simmid} --mutrate ${mr} --max-indel-size 5 --samples ${ss} --migration ${mm}"

      # Add the command to the commands array
      commands+=("$cmd")
    done
  done
done
echo "There are ${#commands[@]} commands to be submitted."

# Submit the commands in batches of 20 using ParaFly
batch_size=1
num_batches=$((${#commands[@]} / batch_size))

echo "Beginning submission of all combinations in parallel batches of ${batch_size}..."
conda activate machina
#for ((i=0; i<$num_batches; i++))
for ((i=0; i<2; i++))
do
  echo "Submitting batch ${i} of ${num_batches}:"
  # Calculate the range of commands to submit in this batch
  start_index=$((i * batch_size))

  # Slice the commands array to get the commands for this batch
  batch_commands=("${commands[@]:${start_index}:${batch_size}}")
  for command in "${batch_commands[@]}"
  do
    echo "${command}" >> "${i}.cmd"
    echo "${command}"
  done
  # Print the ParaFly command to the console
  ParaFly -CPU ${batch_size} -c ${i}.cmd
  rm ${i}.cmd
  rm ${i}.cmd.completed
done

echo "Commands are done, now appending results to one file..."

# Create the main output file with the header
echo "name,mutrate,max_indel_size,num_samples,migration_matrix,true_migrations,inferred_migrations,proportion" > ${parallel_sim_name}_all_output.csv

# Append the second line of each csv file to the main file
for f in sim_results*/comparison*; do tail -n +2 $f >> ${parallel_sim_name}_all_output.csv; done

# Remove all the individual csv files
rm -r sim_results*

echo "Done."

conda deactivate

