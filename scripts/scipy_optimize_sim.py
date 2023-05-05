#!/usr/bin/env python3
import os
import subprocess
from scipy.optimize import minimize_scalar
import pandas as pd
import re
import seaborn as sns
import matplotlib as mpl

#I will set parallel_sim.sh to run 100x simulations for the input mutrate and only return the proportion values from the 
# summary csv as 'proportions = [0.1234 0.2345 0.3213 ...]' format in the last line of terminal output for simplicity with scipy.optimize

def optimize_constant_mutrate(x):
    input_string = f'{x},{x},{x},{x},{x},{x},{x},{x},{x},{x}'
    result = subprocess.run(['bash', 'parallel_sim.sh', '-m', input_string], stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8').split('\n')
    print(f'Input mutrate string : {input_string}')
    proportion_regex = r'(?<=proportions = \[)[\d\.,]+(?=\])'
    match = re.search(proportion_regex, output[-2])
    if match:
        proportions_str = match.group(0)
        proportions = [float(p) if p != '' else 0 for p in proportions_str.split(',')]
        average = sum(proportions)/len(proportions)
        print(f"Average proportion true migrations inferred: {average}")
    else:
        print("No proportion values found in the last line.")
        proportions = []
    results[input_string] = proportions
    print(results)
    return -(average) if proportions else None

def avg_mutrate(mutrate_str):
    mutrates = [float(x) for x in mutrate_str.split(',')]
    return sum(mutrates) / len(mutrates)

#Set output name
output_name="constantMutrate_scipy_results"
output_dir=f'./{output_name}'
results = {}

#Run the scipy optimization
result = minimize_scalar(optimize_constant_mutrate, method='bounded', bounds=[0.05,0.25], options={'maxiter':23})

#Print maximum proportion and corresponding input
max = -result.fun
print(f"Maximum value: {max}")
print(f"Location of maximum: {result.x}")

#Save the result to a file
os.makedirs(output_dir)
with open(f'{output_dir}/result_{output_name}.txt', 'w') as f:
    f.write(f"Status: {result.success}\n")
    f.write(f"Solution: {result.x}\n")
    f.write(f"Function value: {result.fun}\n")
    f.write(f"Number of function evaluations: {result.nfev}\n")

#Save results to csv
results_df = pd.DataFrame.from_dict(results, orient='index')
results_df.reset_index(drop=False, inplace=True)
results_df.rename(columns={'index':'mutrate'}, inplace=True)
results_df['avg_mutrate'] = results_df['mutrate'].apply(avg_mutrate)
results_df.to_csv(f'{output_dir}/data_{output_name}.csv', index=False)

#Make simple avg mutrate vs proportion plot to visualize results
#sns.lineplot(x='avg_mutrate',y=,data=results_df, errorbar=('ci', 95), color='blue')
#plt.tight_layout()
#plt.savefig(f'{output_dir}/plot_{output_name}.png')



