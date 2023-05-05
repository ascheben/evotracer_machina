import sys
import pandas as pd
import numpy as np

# This script takes in the filepath to the simmulated tissues.tsv file, the out machina migration.txt file, 
# the sim name, and the output directory passed from the simulation name
true_migration_filepath = sys.argv[1]
inferred_migration_filepath = sys.argv[2]
sim_name = sys.argv[3]
outdir = sys.argv[4]

true_df = pd.read_csv(true_migration_filepath, sep="\t")
inferred_df = pd.read_csv(inferred_migration_filepath)

### Convert true df format to an interpreted migration format similar to the inferred output by the machina wrapper
tissues = true_df['tissue'].unique()
true_transitions = {}
for from_tissue in tissues:
    for to_tissue in tissues:
        transition_key = f"{from_tissue}:{to_tissue}"
        true_transitions[transition_key] = 0

# loop over each row in the DataFrame
for index, row in true_df.iterrows():
    # get the node and parent_node values
    node = row['node']
    parent_node = row['parent_node']
    
    # if the parent_node value is not null
    if pd.notnull(parent_node):
        # get the tissue labels for the current node and parent node
        node_tissue = row['tissue']
        parent_tissue = true_df.loc[true_df['node'] == parent_node, 'tissue'].iloc[0]
        
        # create a key for the tissue transition and increment the count in the dictionary
        transition_key2 = f"{parent_tissue}:{node_tissue}"
        true_transitions[transition_key2] += 1
       

# convert the dictionary to a DataFrame and sort by the transition counts in descending order
transitions_df = pd.DataFrame(list(true_transitions.items()), columns=['transition', f'{sim_name}_true']).sort_values(f'{sim_name}_true', ascending=False)
transitions_df = transitions_df.set_index('transition')
transitions_df = transitions_df.transpose()


### Extract the inferred migration values form the dataframe and append under the true values in the transitions df above
transitions_df.loc[f'{sim_name}_inferred']=inferred_df.loc[:,inferred_df.columns.str.contains(':')].sum()


# output the comparison DataFrame to a CSV file
transitions_df.reset_index(drop=False, inplace=True)
transitions_df = transitions_df.rename(columns={'index':'name'})
transitions_df.to_csv(f'{outdir}detailed_comparison_inferred_true_migration_{sim_name}.csv', index=False, na_rep='NaN')

