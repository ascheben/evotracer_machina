#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

parallel_sim_name="samplesize2_explore_parameters"
# Set the migration matrix values that will be tested
rare_mm="data/rare_migration_prob_matrix.csv"
equal_mm="data/equal_migration_prob_matrix.csv"
true_mm="data/true_migration_prob_matrix.csv"

# Set the mutation matrix filepaths
#mm_array=(${equal_mm} ${true_mm} ${rare_mm})
mm_array=(${true_mm})

# Set the mutation rates to explore
#mr_array=(0.01 0.025 0.05 0.1 0.2)
mr_array=(0.05)

# Set the max indel sizes to explore
#mi_array=(3 5 10 15 20)
mi_array=(5)

# Set the sample sizes to explore
ss_array=(100 200 300 400 500)
#ss_array=(100)

# Setup headers for recording the input parameters in a csv
mm_header="mutation_matrix"
mr_header="mutation_rate"
mi_header="max_indel_size"
ss_header="sample_size"

# Determine the maximum length of the arrays for input parameter csv writing
max_len=${#mm_array[@]}
if [ ${#mr_array[@]} -gt $max_len ]; then
  max_len=${#mr_array[@]}
fi
if [ ${#mi_array[@]} -gt $max_len ]; then
  max_len=${#mi_array[@]}
fi
if [ ${#ss_array[@]} -gt $max_len ]; then
  max_len=${#ss_array[@]}
fi

# Write the array values to a CSV file
echo "${mr_header},${mi_header},${ss_header},${mm_header}" > input_ranges_${parallel_sim_name}.csv

for ((i=0; i<$max_len; i++)); 
  do
    mr=${mr_array[$i]:-}
    mi=${mi_array[$i]:-}
    ss=${ss_array[$i]:-}
    mm=${mm_array[$i]:-}
    echo "${mr},${mi},${ss},${mm}" >> input_ranges_${parallel_sim_name}.csv
  done

# Create an array to store all the commands
commands=()

echo "Setting up all parameter combinations..."
# Loop over all combinations of the three arrays and add the command to the commands array
for mm in ${mm_array[@]}
do
  for mr in ${mr_array[@]}
  do
    for mi in ${mi_array[@]}
    do
      for ss in ${ss_array[@]}
      do
        # Construct the command to submit with ParaFly
        mm_type=$(echo $mm | awk -F/ '{print $2}' | awk -F_ '{print $1}')
        mr_value=$(echo $mr | awk -F. '{print $2}')
        simmid="sim_${mr_value}_${mi}_${ss}_${mm_type}"
        cmd="./sim_full_pipeline.sh --out ${simmid} --mutrate ${mr} --max-indel-size ${mi} --samples ${ss} --migration ${mm}"
        # Add the command to the commands array
        commands+=("$cmd")
      done
    done
  done
done
echo "There are ${#commands[@]} commands to be submitted."

# Create the main output file with the header
echo "name,mutrate,max_indel_size,num_samples,migration_matrix,true_migrations,inferred_migrations,proportion" > output_all_${parallel_sim_name}.csv

# Submit the commands in batches of 20 using ParaFly
batch_size=5
batches=$((${#commands[@]} / batch_size))
num_batches=$(echo "scale=0; ($batches+0.5)/1" | bc)

echo "Beginning submission of all combinations in parallel batches of ${batch_size}..."
conda activate machina
for ((i=0; i<$num_batches; i++))
do
  echo "Submitting batch $((i+1))/$((num_batches+1)):"
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
  # Append the second line of each csv file to the main file
  for f in sim_results*/comparison*; do tail -n +2 $f >> output_all_${parallel_sim_name}.csv; done
  # Remove all the individual csv files
  rm -r sim_results*
done

echo "All batches are done."

conda deactivate

par_results="${parallel_sim_name}_parallel_sim_results"
mkdir ${par_results}
mv *${parallel_sim_name}.csv ${par_results}/

