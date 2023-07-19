import sys
import pandas as pd
import numpy as np
import os

# This script takes in the filepath to the simmulated tissues.tsv file, the out machina migration.txt file, 
# the sim name, and the output directory passed from the simulation name
true_migration_filepath = sys.argv[1]
inferred_migration_filepath = sys.argv[2]
sim_name = sys.argv[3]
outdir = sys.argv[4]
migration_matrix_filepath = sys.argv[5]

true_df = pd.read_csv(true_migration_filepath, sep="\t")
if os.path.exists(inferred_migration_filepath):
    inferred_df = pd.read_csv(inferred_migration_filepath)

migration_matrix = pd.read_csv(migration_matrix_filepath)

### Convert true df format to an interpreted migration format similar to the inferred output by the machina wrapper
#tissues = true_df['tissue'].unique()
tissues = migration_matrix.columns[1:].values
true_transitions = {}
migration_keys = []
for from_tissue in tissues:
    for to_tissue in tissues:
        transition_key = f"{from_tissue}:{to_tissue}"
        true_transitions[transition_key] = 0
        if (from_tissue != to_tissue) and (f"{from_tissue}:{to_tissue}" not in migration_keys):
            migration_keys.append(f"{from_tissue}:{to_tissue}")

# loop over each row in the DataFrame
for index, row in true_df.iterrows():
    # get the node and parent_node values
    node = row['node']
    parent_node = row['parent_node']
    
    # if the parent_node value is not null
    if pd.notnull(parent_node):
        # get the tissue labels for the current node and parent node
        node_tissue = row['tissue']
        #parent_tissue = true_df.loc[true_df['node'] == parent_node, 'tissue'].iloc[0]      ### referencing parent nodes causes a problem with downsampled tree tissues tsv output. I instead added a parent tissue column to tissue output
        parent_tissue = row['parent_tissue']
        
        # create a key for the tissue transition and increment the count in the dictionary
        transition_key2 = f"{parent_tissue}:{node_tissue}"
        true_transitions[transition_key2] += 1
       

# convert the dictionary to a DataFrame and sort by the transition counts in descending order
transitions_df = pd.DataFrame(list(true_transitions.items()), columns=['transition', f'{sim_name}_true']).sort_values(f'{sim_name}_true', ascending=False)
transitions_df = transitions_df.set_index('transition')
transitions_df = transitions_df.transpose()


### Extract the inferred migration values form the dataframe and append under the true values in the transitions df above
if os.path.exists(inferred_migration_filepath):
    transitions_df.loc[f'{sim_name}_inferred']=inferred_df.loc[:,inferred_df.columns.str.contains(':')].sum()
    transitions_df['migration_matrix'] = [migration_matrix_filepath, migration_matrix_filepath]       # add migration matrix value for identification when running in parallel
else:
    transitions_df['migration_matrix'] = [migration_matrix_filepath]
transitions_df = transitions_df.fillna(0)

# calculate specific metastatic path aware true positives, false positives, and false negative counts with statistics for inferred vs. true data
if os.path.exists(inferred_migration_filepath):
    migration_df = transitions_df[migration_keys]
    tp = 0
    fp = 0
    fn = 0
    for column in migration_df.columns:
        if migration_df.loc[f'{sim_name}_true', column] >=  migration_df.loc[f'{sim_name}_inferred', column]:
            tp += migration_df.loc[f'{sim_name}_inferred', column]
            fn += (migration_df.loc[f'{sim_name}_true', column] - migration_df.loc[f'{sim_name}_inferred', column])
        else:
            tp += migration_df.loc[f'{sim_name}_true', column]
            fp += (migration_df.loc[f'{sim_name}_inferred', column] - migration_df.loc[f'{sim_name}_true', column])
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * (precision * recall) / (precision + recall)
    stats_data = {
        'machina_sim_name': [sim_name],
        'migration_matrix': [migration_matrix_filepath],
        'true_positives': [tp],
        'false_positives': [fp],
        'false_negatives': [fn],
        'precision': [precision],
        'recall': [recall],
        'f1_score': [f1_score]
    }

    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv(f'{outdir}statistics_pathAwareComparison_inferred_true_migration_{sim_name}.csv', index=False)

# output the comparison DataFrame to a CSV file
transitions_df.reset_index(drop=False, inplace=True)
transitions_df = transitions_df.rename(columns={'index':'name'})
transitions_df.to_csv(f'{outdir}detailed_comparison_inferred_true_migration_{sim_name}.csv', index=False, na_rep=0)
    

