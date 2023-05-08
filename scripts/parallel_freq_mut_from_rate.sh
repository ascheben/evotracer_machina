#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

#mutrate_array=(0.01 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.250 0.275 0.3 0.325 0.350 0.375 0.4 0.425 0.45 0.475 0.5)
mutrate_array=(0.525 0.550 0.575 0.6 0.625 0.650 0.675 0.7 0.725 0.750 0.775 0.8 0.825 0.850 0.875 0.9 0.925 0.950 0.975 1)
output_name="ext12_experimental_mutrate_search"

conda activate simulate

j=2100
for mr in ${mutrate_array[@]};
do
    for ((rep=0; rep<100; rep++))   ### Can set the amount of repeat simulations with the same parameters
      do
      output="run${j}"
      cmd="python scripts/freq_mut_from_rate.py ${mr} ${output}"
      echo "${cmd}" >> "${output_name}.cmd"
      j=$((j+1))
      done
done

batch_size=20

ParaFly -CPU ${batch_size} -c ${output_name}.cmd
  rm ${output_name}.cmd
  rm ${output_name}.cmd.completed

echo "output_name,mut_rate,num_mut_leaves,sample_size,freq_mut_outcome" >> ${output_name}_freq_mut_from_rate.csv
cat *_freq_mut_from_rate_data.csv >> ${output_name}_freq_mut_from_rate.csv
rm *_freq_mut_from_rate_data.csv

conda deactivate


