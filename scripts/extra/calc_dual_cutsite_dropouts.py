import pandas as pd
import sys

### This script takes in the simulator output of indel character matrix tsv, mutations tsv, the simmid name, and a desired output directory
### It then calculates the number of samples/rows in the indel character matrix that have two or more mutations with the same age, indicating two or more cutsites being cut at the same time
indel_matrix_filepath=sys.argv[1]
mutation_key_filepath=sys.argv[2]
simmid=sys.argv[3]
outdir=sys.argv[4]

# Read in the two TSV files as pandas dataframes
df1 = pd.read_csv(indel_matrix_filepath, sep='\t', index_col=0)
df2 = pd.read_csv(mutation_key_filepath, sep='\t', header=None)

# Create a dictionary to map the keys to their corresponding ages
key_to_age = dict(zip(df2[0], df2[2]))

# Create a new dataframe to hold the results
results_df = pd.DataFrame(columns=['sample', 'total_mutations', 'num_duplicate_ages', 'proportion_duplicate_cuts'])

# Loop through each row in the first dataframe
for sample, row in df1.iterrows():
    # Create a list to hold the ages of the non-zero keys in this row
    ages = []
    
    # Loop through each non-zero value in the row
    for value in row:
        if value != 0:
            # If the value is non-zero, add the corresponding age to the set
            ages.append(key_to_age[value])
    
    # Count the number of duplicate ages in the set
    total_mutations=len(ages)
    num_duplicate_ages = len(ages) - len([age for age in ages if ages.count(age) == 1])
    proportion_dup_cuts=round((num_duplicate_ages/len(ages)), 4)
    # Append the results for this sample to the results dataframe
    results_df = results_df.append({'sample': sample, 'total_mutations': total_mutations, 'num_duplicate_ages': num_duplicate_ages, 'proportion_duplicate_cuts': proportion_dup_cuts}, ignore_index=True)

# Print the results dataframe
results_df['sample'] = results_df['sample'].astype(int)
results_df['total_mutations'] = results_df['total_mutations'].astype(int)
results_df['num_duplicate_ages'] = results_df['num_duplicate_ages'].astype(int)
results_df.to_csv(f"{outdir}{simmid}_dual_cutsite_dropout_results.csv", index=False)
