import sys
import pandas as pd
import numpy as np
import re

migration_filepath = sys.argv[1]
outprefix = sys.argv[2]
desired_primary_tissue = sys.argv[3]

#migration_filepath = "test_mig_cp_output/test_mig_migration.txt"
#outprefix = "data/true"
#desired_primary_tissue = 'PRL'


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

# Normalize probabilities for each prior tissue to sum to 1 for downstream random choice function
for prior_tissue, transition_probabilities in prob_dict.items():
    probabilities = list(transition_probabilities.values())
    normalized_probabilities = np.round(probabilities, 2) / np.sum(np.round(probabilities, 2))
    new_transition_probabilities = dict(zip(transition_probabilities.keys(), np.round(normalized_probabilities,2)))
    prob_dict[prior_tissue] = new_transition_probabilities
    
# Re-order the probability dictionary to have the desired primary tissue as the first key
rest_keys = set(prob_dict.keys()) - {desired_primary_tissue}
new_prob_dict = {desired_primary_tissue: prob_dict[desired_primary_tissue]}
for key in rest_keys:
    new_prob_dict[key] = prob_dict[key]

# Write the transition matrix as an output csv file
output_df = pd.DataFrame.from_dict(new_prob_dict, orient="index")
output_df = output_df[output_df.index]
output_df.to_csv(outprefix + '_migration_prob_matrix.csv')


    



