import pandas as pd
import test_config
import argparse
import difflib
import json
import os

def convert_to_percentage(a, b):
    return round((a/b)*100, 2)

def diff_letters(a,b):
    return len([li for li in difflib.ndiff(a, b) if li[0] != ' '])

def output_mutation_dicts_list(file_path, mut_dicts):
    with open(file_path, 'w') as file:
        for mut_dict in mut_dicts:
            json.dump(mut_dict, file)
            file.write('\n')   

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--output-dir", default="test_output")
parser.add_argument("-a", "--output-all-mutations", action="store_true", help="Output all simulated and evotracer ASV cut-site mutation mappings to file")
parser.add_argument("-s", "--output-unidentified-sim-mutations", action="store_true", help="Output all simulated ASV mutations that are not identified by EvoTraceR")
parser.add_argument("-i", "--output-identified-sim-mutations", action="store_true", help="Output all simulated ASV mutations that are identified by EvoTraceR")
parser.add_argument("-m", "--map-incorrect-evo-mutations", action="store_true", help="EvoTraceR ASV mutations without an exact match in the simulated data are mapped to the most similar ASV mutations in the simulated data. Results are output to file.")
parser.add_argument("-p", "--pos-error", type=int, default=10, help="Number of position shifts allowed in mutations when finding similar ASVs")
parser.add_argument("-c", "--indel-char-error", type=int, default=1, help="Number of character mismatches allowed in mutations when finding similar ASVs")
args = vars(parser.parse_args())

output_dir_path = f'{os.getcwd()}/{args["output_dir"]}'
if not os.path.exists(output_dir_path):
    os.mkdir(output_dir_path)

# prepare data
sim_mutation_matrix = pd.read_csv(test_config.indel_character_matrix_file_path, sep='\t')
sim_mutation_matrix = sim_mutation_matrix.iloc[:, 1:]
cut_positions = pd.read_csv(test_config.cut_positions_file_path, sep='\t', header=None)
cut_positions = {c[0]: c[1] for c in cut_positions.values}
mutations = pd.read_csv(test_config.mutations_file_path, sep='\t', header=None)
mutations = {m[0]: m[1] for m in mutations.values}
sim_mutation_matrix = sim_mutation_matrix.rename(cut_positions, axis=1)
sim_mutation_matrix = sim_mutation_matrix.drop_duplicates()
evo_df = pd.read_csv(test_config.asv_stat_file_path)
evo_df = evo_df.loc[evo_df['asv_names'] != test_config.ref_name] # remove ref sequence
evo_asv_seqs = evo_df['seq'].unique()
ref_seq = test_config.ref_seq

# get dictionaries of cut site mutation mappings for each unique simulated sequence
sim_cut_site_mutation_mappings = []
for index, row in sim_mutation_matrix.iterrows():
    row = row.to_frame()
    row = row[(row != 0).all(1)] # remove non-mutated entries
    seq_mutations = row.to_dict()[index]
    for pos,mut_id in seq_mutations.items():
        seq_mutations[pos] = mutations[mut_id]
    sim_cut_site_mutation_mappings.append(seq_mutations)

if args["output_all_mutations"]: 
    output_mutation_dicts_list(f'{output_dir_path}/simulated_cut_site_mutation_mappings.txt', sim_cut_site_mutation_mappings)
            
# generate simulated ASVs using the cut site mutation mappings and the ref sequence
sim_asv_seqs = []
total_simulated_mutations = 0
for simulated_asv in sim_cut_site_mutation_mappings:
    total_simulated_mutations += len(simulated_asv)
    mutated_seq = ref_seq
    shift = 0
    for pos, alt_seq in simulated_asv.items():
        pos += shift - 1
        if '-' in alt_seq:
            mutated_seq = mutated_seq[:pos] + mutated_seq[pos + len(alt_seq):]
            shift -= len(alt_seq)
        else:
            mutated_seq = mutated_seq[:pos] + alt_seq + mutated_seq[pos:]
            shift += len(alt_seq)
    sim_asv_seqs.append(mutated_seq)

# # Test Case #1: Percentage of Simulated ASVs identified
matching_asvs = set(sim_asv_seqs) & set(evo_asv_seqs)
print("Compare simulated sequences to asv_stat seq column values:")
print(f'Simulated data has {len(sim_asv_seqs)} ASVs')
print(f'EvoTracer found {len(evo_asv_seqs)} ASVs')
print(f'Percentage of Simulated ASVs found: {convert_to_percentage(len(matching_asvs), len(sim_asv_seqs))}%\n')

# # Test Case #2: Percentage of mutations identified, calculated using mutation, 
# # cut site, and character indel matrix data
evo_df = evo_df.drop_duplicates(subset=['asv_names', 'mutation_type', 'start', 'alt_seq'])

# get evotracer mutations from asv_stat data
evo_asv_muts_dict = {}
for index, row in evo_df.iterrows():
    if row.asv_names not in evo_asv_muts_dict:
        evo_asv_muts_dict[row.asv_names] = {}
    evo_asv_muts_dict[row.asv_names][int(row.start)] = row.alt_seq

if args["output_all_mutations"]:
    with open(f'{output_dir_path}/evotracer_cut_site_mutation_mappings.txt', 'w') as file:
        for asv, mutations in evo_asv_muts_dict.items():
            file.write(f'{asv}: {json.dumps(mutations)}\n')

# find simulated set of mutations matching evotracer set of mutations, otherwise find most similar mutations
total_evotracer_mutations = 0
total_identified_mutations = 0
identified_mutations = []
pos_shift_err = args["pos_error"]
seq_char_err = args["indel_char_error"]
if args["map_incorrect_evo_mutations"]: open(f'{output_dir_path}/incorrect_evotracer_mutations_mappings.txt', "w")
for asv, evo_mutations in evo_asv_muts_dict.items():
    evo_mutations = dict(sorted(evo_mutations.items())) # sort mutations in dictionary by position
    total_evotracer_mutations += len(evo_mutations)
    # identify matching mutations
    if evo_mutations in sim_cut_site_mutation_mappings:
        total_identified_mutations += len(evo_mutations)
        identified_mutations.append(evo_mutations)
    # identify most similar set of mutations in simulated data compared to evotracer generated ASV mutations
    elif args["map_incorrect_evo_mutations"]:
        # filter simulated mutation set with same number of mutations as evotracer ASV
        # similar_mutations_list = [x for x in sim_cut_site_mutation_mappings if len(x) == len(evo_mutations)]
        similar_mutations = []
        for simulated_mutations in sim_cut_site_mutation_mappings:
            simulated_mutations = dict(sorted(simulated_mutations.items()))
            # filter simulated mutation set by positions and indel characters
            if len(simulated_mutations) == len(evo_mutations):
                is_similar = True
                for (simulated_pos, simulated_alt_seq), (evo_pos, evo_alt_seq) in zip(simulated_mutations.items(), evo_mutations.items()):
                    if simulated_pos < evo_pos - pos_shift_err or simulated_pos > evo_pos + pos_shift_err or diff_letters(simulated_alt_seq, evo_alt_seq) > seq_char_err:
                        is_similar = False
                if is_similar:
                    similar_mutations.append(simulated_mutations)
        incorrect_evo_mut_file_path = f'{output_dir_path}/incorrect_evotracer_mutations_mappings.txt'
        with open(incorrect_evo_mut_file_path, "a") as file:
            file.write(f"EvoTraceR {asv}: {evo_mutations}\n")
            file.write("Most similar simulated ASVs:\n")
            for mut_dict in similar_mutations:
                json.dump(mut_dict, file)
                file.write('\n')
            file.write('\n')   
                    
print("Compare simulated sequence mutations to asv_stat start and alt_seq column values:")
print(f"Percentage of mutations found: {convert_to_percentage(total_evotracer_mutations, total_simulated_mutations) if total_evotracer_mutations <= total_simulated_mutations else 100.0}%")
print(f"Percentage of mutations with accurate positions and indel characters: {convert_to_percentage(total_identified_mutations, total_simulated_mutations)}%")

unidentified_sim_mutations = [dict for dict in sim_cut_site_mutation_mappings if dict not in evo_asv_muts_dict.values()]

if args["output_unidentified_sim_mutations"]:
    output_mutation_dicts_list(f'{output_dir_path}/unidentified_simulated_mutations.txt', unidentified_sim_mutations)

if args["output_identified_sim_mutations"]:
    output_mutation_dicts_list(f'{output_dir_path}/identified_simulated_mutations.txt', identified_mutations)
