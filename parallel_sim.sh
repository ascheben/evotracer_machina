#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

parallel_sim_name="trueMM_subpopulationMigrations"

# Set the mutation rates to explore
# mr1=(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
# mr2=(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0)
# mr3=(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0,0)
# mr4=(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0,0,0)
# mr5=(0.1,0.1,0.1,0.1,0.1,0.1,0,0,0,0)
# mr6=(0.1,0.1,0.1,0.1,0.1,0,0,0,0,0)
# mr7=(0.1,0.1,0.1,0.1,0,0,0,0,0,0)
# mr8=(0.1,0.1,0.1,0,0,0,0,0,0,0)
# mr9=(0.1,0.1,0,0,0,0,0,0,0,0)
# mr10=(0.1,0,0,0,0,0,0,0,0,0)
# mr11=(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
# mr12=(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
# mr13=(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
# mr14=(0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15)
# mr15=(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)
# mr16=(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3)
# mr17=(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
# mr18=(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7)
# mr19=(1,1,1,1,1,1,1,1,1,1)
# mr20=(1,0.01,1,0.01,1,0.01,1,0.01,1,0.01)
# mr21=(1,1,1,1,0.1,0.1,0.1,0.01,0.01,0.01)
# mr22=(1,1,1,1,1,0.01,0.01,0.01,0.01,0.01)
# mr23=(1,1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
# mr24=(1,0.1,1,0.1,1,0.1,1,0.1,1,0.1)
# mr25=(1,0.01,0.01,1,0.01,0.01,1,0.01,0.01,1)
# mr26=(1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
# mr27=(0.1,0.01,0.1,0.01,0.1,0.01,0.1,0.01,0.1,0.01)
# mr_array=("$mr1" "$mr2" "$mr3" "$mr4" "$mr5" "$mr6" "$mr7" "$mr8" "$mr9" "$mr10" "$mr11" "$mr12" "$mr13" "$mr14" "$mr15" "$mr16" "$mr17" "$mr18" "$mr19" "$mr20" "$mr21" "$mr22" "$mr23" "$mr24" "$mr25" "$mr26" "$mr27")

### mr input for mied strategies with an average of 0.1264
mr1=(0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264,0.1264)
mr2=(0.3173305,0.0009823,0.00000008,0.00000162,0.00000335,0.18609872,0.00000254,0,0,0)
mr_array=("$mr1")

### Use below to take input mutrate for scipy.optimize script
#while [[ "$#" -gt 0 ]]; do
#    case $1 in
#        -m|--mutrate) mr_array="$2"; shift ;;
#    *) echo "Unknown parameter passed: $1"; echo "Usage: parallel_sim.sh -m <10 comma seperated values 0 to 1>" ; exit 1 ;;
#    esac
#    shift
#done ### remove to above if not doing scipy for mutrate. Input mutrate can be manual for manual simulations not with scipy.

# Set the sample sizes to explore
#ss_array=(100 250 500 750 1000)
ss_array=(100)

# Set the migration matrix values that will be tested
ultra_rare_mm="data/ultrarare_migration_prob_matrix.csv"
rare_mm="data/rare_migration_prob_matrix.csv"
moderate_mm="data/moderate_migration_prob_matrix.csv"
true_mm="data/true_migration_prob_matrix.csv"
high_mm="data/high_migration_prob_matrix.csv"
equal_mm="data/equal_migration_prob_matrix.csv"
#mm_array=(${rare_mm} ${equal_mm} ${moderate_mm} ${high_mm} ${true_mm})
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
      for ((rep=0; rep<100; rep++))   ### Can set the amount of repeat simulations with the same parameters
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
echo "name,mutrate,num_samples,migration_matrix,num_uniq_mutations,uniq_mut_per_site,total_mut_sites,avg_mut_sites_per_sample,avg_proportion_mut_sites,total_dropout,dropout_per_sample,avg_mutation_age,true_migrations,inferred_migrations,proportion" > output_all_${parallel_sim_name}.csv

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


### Use below to print all proportions data and remove sim data for use with scipy optimize
### Remove this if not using scipy since this deletes all output and prints to the terminal only what scipy needs
#proportions_scipy=$(tail -n +2 ${par_results}/output_all_${parallel_sim_name}.csv | cut -d ',' -f 15 | sed ':a;N;$!ba;s/\n/,/g')
#echo "proportions = [$proportions_scipy]"
#rm -r ${par_results}/


