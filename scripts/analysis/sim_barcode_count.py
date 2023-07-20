#import pandas as pd
import sys

""" indel_character_matrix_filepath = sys.argv[1]
mutation_key_filepath = sys.argv[2]

# Read the first matrix file into a DataFrame
barcode_df = pd.read_csv(indel_character_matrix_filepath, sep='\t', index_col=0)
mutations_df = pd.read_csv(mutation_key_filepath, sep='\t', index_col=False, usecols=[0,1], header=None, names=['Index', 'String'])
mutations_df = mutations_df.astype(str)
mutations_df = mutations_df.astype(str)
mutations_dict = mutations_df[['Index', 'String']].set_index('Index').to_dict(orient='dict')['String']
mutations_dict['-1'] = ''
mut_barcode_df = barcode_df.replace(mutations_dict)
mut_barcode_df = mut_barcode_df.astype(str)

concatenated_mutations = mut_barcode_df.apply(lambda row: ''.join(row), axis=1)
unique_strings_count = len(concatenated_mutations.unique())
print(unique_strings_count)
 """

fasta_filepath = sys.argv[1]

# Open the FASTA file
with open(fasta_filepath, 'r') as file:
    # Initialize a set to store unique character rows
    unique_rows = set()

    # Iterate over each line in the file
    for line in file:
        # Remove leading/trailing whitespaces
        line = line.strip()

        # Skip the reference lines starting with ">"
        if line.startswith(">"):
            continue

        # Add the character row to the set
        unique_rows.add(line)

# Count the number of unique character rows
num_unique_rows = len(unique_rows)

# Print the result
print(num_unique_rows)



