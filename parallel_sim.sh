#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

parallel_sim_name="test_explore_parameters"

# Set the mutation rates to explore
mr1=(0.1,0.1,0.1,0.1,0.1,0,0,0,0,0)
mr2=(0.1,0.1,0.1,0.2,0.2,0,0,0,0,0)
mr3=(0.2,0.2,0.2,0.3,0.3,0,0,0,0,0)
#mr_array=("$mr1" "$mr2" "$mr3")
mr_array=(0.1,0.1,0.1,0.1,0.1,0,0,0,0,0)

# Set the sample sizes to explore
#ss_array=(100 250 500 750 1000)
ss_array=(100)

# Set the migration matrix values that will be tested
rare_mm="data/rare_migration_prob_matrix.csv"
equal_mm="data/equal_migration_prob_matrix.csv"
true_mm="data/true_migration_prob_matrix.csv"
#mm_array=(${equal_mm} ${true_mm} ${rare_mm})
mm_array=(${true_mm})

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
      for ((rep=0; rep<22; rep++))   ### Can set the amount of repeat simulations with the same parameters
      do
      simmid="sim${j}"
      cmd="./sim_full_pipeline.sh --out ${simmid} --mutrate ${mr} --samples ${ss} --migration ${mm}"
      commands+=("$cmd")
      j=$((j+1))
      done
    done
  done
done
echo "There are ${#commands[@]} commands to be submitted."

# Create the main output file with the header
echo "name,mutrate,num_samples,migration_matrix,num_mutations,true_migrations,inferred_migrations,proportion" > output_all_${parallel_sim_name}.csv


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

