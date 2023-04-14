import sys
import pandas as pd
import numpy as np
import re

migration_filepath = sys.argv[1]
outprefix = sys.argv[2]

#migration_filepath = "test_mig_cp_output/test_mig_migration.txt"
#outprefix = "true"


# Read in migration file and extract tissue names
migration_df = pd.read_csv(migration_filepath)
tissue_cols = [col for col in migration_df.columns if re.match(r'^H_', col)]
tissues = [col.split('_', maxsplit=1)[1] for col in tissue_cols]

# Obtain dictionaries of transitions and totals to calculate probabilities from a known prior tissue
transition_dict = {}
total_dict = {}

for primary_tissue in tissues:
    transition_cols = [col for col in migration_df.columns if col.startswith(primary_tissue + ':')]
    total_tissue = 0
    for transition in transition_cols:
        transition_dict[transition] = np.sum(migration_df[transition].values)
        total_tissue = total_tissue + np.sum(migration_df[transition].values)
    total_dict[primary_tissue] = total_tissue

# Calculate the transition matrix as a probability dictionary
prob_dict = {}

#for key in transition_dict:
#    prior_tissue, new_tissue = key.split(':')
#    prob_dict[key] = transition_dict[key] / total_dict[prior_tissue]

for key in transition_dict:
    prior_tissue, new_tissue = key.split(':')
    if prior_tissue in total_dict:
        if prior_tissue not in prob_dict:
            prob_dict[prior_tissue] = {new_tissue: transition_dict[key] / total_dict[prior_tissue]}
        else:
            prob_dict[prior_tissue][new_tissue] = transition_dict[key] / total_dict[prior_tissue]
    
# Write the transition matrix as an output csv file
output_df = pd.DataFrame.from_dict(prob_dict, orient="index")
output_df.to_csv(outprefix + '_migration_prob_matrix.csv')


    



