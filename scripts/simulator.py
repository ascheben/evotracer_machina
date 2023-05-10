#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
import math
import random
from ete3 import Tree
import cassiopeia as cas
from cassiopeia.data import CassiopeiaTree

# monkey patch Cassiopeia to output mutation age
def overlay_data_drop(self, tree: CassiopeiaTree):
    """Overlays Cas9-based lineage tracing data onto the CassiopeiaTree.
    Args:
        tree: Input CassiopeiaTree
    """
    if self.random_seed is not None:
        np.random.seed(self.random_seed)

    # create state priors if they don't exist.
    # This will set the instance's variable for mutation priors and will
    # use this for all future simulations.
    if self.mutation_priors is None:
        self.mutation_priors = {}
        probabilites = [
            self.state_generating_distribution()
            for _ in range(self.number_of_states)
        ]
        Z = np.sum(probabilites)
        for i in range(self.number_of_states):
            self.mutation_priors[i + 1] = probabilites[i] / Z

    number_of_characters = self.number_of_cassettes * self.size_of_cassette

    # initialize character states
    character_matrix = {}
    for node in tree.nodes:
        character_matrix[node] = [-1] * number_of_characters

    for node in tree.depth_first_traverse_nodes(tree.root, postorder=False):

        if tree.is_root(node):
            character_matrix[node] = [0] * number_of_characters
            continue

        parent = tree.parent(node)
        life_time = tree.get_time(node) - tree.get_time(parent)

        character_array = character_matrix[parent]
        open_sites = [
            c
            for c in range(len(character_array))
            if character_array[c] == 0
        ]

        new_cuts = []
        for site in open_sites:
            mutation_rate = self.mutation_rate_per_character[site]
            mutation_probability = 1 - (np.exp(-life_time * mutation_rate))

            if np.random.uniform() < mutation_probability:
                new_cuts.append(site)

        # collapse cuts that are on the same cassette
        cuts_remaining = new_cuts
        if self.collapse_sites_on_cassette and self.size_of_cassette > 1:
            character_array, cuts_remaining = self.collapse_sites(
                character_array, new_cuts
            )
        # introduce new states at cut sites

        character_array = self.introduce_states(
            character_array, cuts_remaining
        )
        if len(cuts_remaining) >0:
            timestamp = round(tree.get_time(node),6)
            for i in cuts_remaining:
                character_array[i] = str(character_array[i]) + "_" + str(timestamp)
            #print(node,tree.get_time(node),character_array, cuts_remaining)
        # silence cassettes
        silencing_probability = 1 - (
            np.exp(-life_time * self.heritable_silencing_rate)
        )
        character_array = self.silence_cassettes(
            character_array,
            silencing_probability,
            self.heritable_missing_data_state,
        )

        character_matrix[node] = character_array

    # apply stochastic silencing
    for leaf in tree.leaves:
        character_matrix[leaf] = self.silence_cassettes(
            character_matrix[leaf],
            self.stochastic_silencing_rate,
            self.stochastic_missing_data_state,
        )

    tree.set_all_character_states(character_matrix)

def mutate_seq(seq,pos,mut_bases):
    '''
    Take sequence and position and randomly introduce an
    indel at the given position. Returns the mutated sequence.
    '''
    mut_size = len(mut_bases)
    if '-' not in mut_bases:
        mut_seq = seq[:pos] + mut_bases + seq[pos:]
    else:
        # pos is 1-based seq idx
        mut_seq = seq[:pos] + mut_bases + seq[pos+mut_size:]
        #mut_seq = seq[:pos-1] + mut_bases + seq[pos-1+mut_size:]
    return mut_seq

def generate_indel():
    '''
    Take sequence and position and randomly introduce an
    indel at the given position. Returns the mutated sequence.
    '''
    # hard limit at 100bp size
    max_size = 100
    del_prob = 0.5
    ins_prob = 0.5
    #max_size = 15
    #max_size = 3
    mut_size = math.ceil(np.random.exponential(6,1)[0])
    while mut_size > max_size:
        mut_size = math.ceil(np.random.exponential(6,1)[0])
    mut_type = random.choices(population=['ins','mut'],weights=[ins_prob,del_prob],k=1)[0]
    #mut_size = random.randint(1, 10)
    if mut_type == 'ins':
        #generate random kmer of mut_size bp
        mut_bases = ''.join(random.choices(['A','C','T','G'], k=mut_size))
    else:
        mut_bases = mut_size * "-"
    return mut_bases

def sim_chars(tree,mut_rate,num_cuts):
    np.random.seed(seed=None)
    lt_sim = cas.sim.Cas9LineageTracingDataSimulator(
        number_of_cassettes = num_cuts,
        size_of_cassette = 1,
        mutation_rate = mut_rate, # can also be a list of rates of equal len to num_cuts
        state_generating_distribution = lambda: np.random.exponential(1e-5),
        number_of_states = 100,
        state_priors = None,
        heritable_silencing_rate = 0, #heritable_silencing_rate = 9e-4,
        stochastic_silencing_rate = 0, #stochastic_silencing_rate = 0.1,
        heritable_missing_data_state = -1,
        stochastic_missing_data_state = -1,
    )
    lt_sim.overlay_data(tree)
    character_matrix = tree.character_matrix
    #character_matrix = ground_truth_tree.character_matrix
    return character_matrix

def assign_tissue_labels(cas_tree,trans_mat):
    tree = cas_tree
    tree_labeled = tree.copy()
    tissue_array = np.array(list(trans_mat.keys()))
    tis_labels = {}
    parent_nodes = {}
    parent_tissues = {}
    migration_labels = {}
    migration_scaling_boundary = 2.5
    migration_scaling_factor = 0.1
    
    # Set the root node label to the first tissue in the list
    root = tree.root
    tis_labels[root] = tissue_array[0]
    
    # Traverse the tree and assign tissue labels to each node
    for node in tree.nodes:
        # Skip the root node since it has already been labeled
        if node == root:
            continue

        # Determine the probability of changing tissue label
        prev_tissue = tis_labels[tree.parent(node)]
        prob_vec = np.array([trans_mat[prev_tissue][t] for t in tissue_array])

        # Make tissue label decision of the current node and assign the label
        if tree.is_leaf(node) == False:
            # Scale the migration probability to be less at earlier points in the tree
            age = tree.get_time(node)
            if age < migration_scaling_boundary:
                tissue_index = np.where(tissue_array == prev_tissue)
                prob_vec = prob_vec * migration_scaling_factor
                sum_except_index = np.sum(np.delete(prob_vec, tissue_index))
                prob_vec[tissue_index] = 1 - sum_except_index
                tis_labels[node] = np.random.choice(tissue_array, p=(prob_vec))
            else:
                tis_labels[node] = np.random.choice(tissue_array, p=prob_vec)
        # Label leaves not by probability but by their parent
        if tree.is_leaf(node) == True:
            tis_labels[node] = prev_tissue
        
        parent_nodes[node] = tree.parent(node)
        parent_tissues[node] = prev_tissue
        
        if tis_labels[node] == prev_tissue:
            migration_labels[node] = False
        if tis_labels[node] != prev_tissue:
            migration_labels[node] = True
    
    # Make a concatenated key and value dictionary to relabel nodes in copied tree
    node_tis_map = {}
    for key, value in tis_labels.items():
        new_value = key + "_" + value
        node_tis_map[key] = new_value
    
    tree_labeled.relabel_nodes(node_tis_map)

    # Add tissue labels to pandas df to add column labeling leaves
    tissues_df = pd.DataFrame.from_dict(tis_labels, orient='index', columns=['tissue'])
    tissues_df.index.names = ['node']
    tissues_df.reset_index(drop=False, inplace=True)
    tissues_df['leaves'] = [tree.is_leaf(x) for x in tissues_df['node'].values]
    tissues_df['parent_node'] = tissues_df['node'].map(parent_nodes)
    tissues_df['parent_tissue'] = tissues_df['node'].map(parent_tissues)
    tissues_df['migration_event'] = tissues_df['node'].map(migration_labels)

    # Change Casseiopeia tree to ETE tree to retain internal node labels when writing newick later
    connections = tree_labeled.edges
    tree_labeled = Tree.from_parent_child_table(connections)
    
    return tissues_df, tree_labeled

def subsample_tissue_labels(tissue_df,subsampled_tree):
    tissue_df = tissue_df.copy()
    sub_tree = subsampled_tree.copy()

    # re-label the subsampled tree with the tissue label concatenated
    subsampled_nodes = sub_tree.nodes
    tissue_df_subsampled = tissue_df[tissue_df['node'].isin(subsampled_nodes)]
    node_tissue_dict = {key: f'{key}_{value}' for key, value in zip(tissue_df_subsampled['node'], tissue_df_subsampled['tissue'])}
    sub_tree.relabel_nodes(node_tissue_dict)

    # Change Casseiopeia tree to ETE tree to retain internal node labels when writing newick later
    connections = sub_tree.edges
    sub_tree = Tree.from_parent_child_table(connections)

    return tissue_df_subsampled, sub_tree


## PARAMETERS ##
outprefix = sys.argv[1]
m = sys.argv[2]
sample_num = int(sys.argv[3])
if sample_num > 10000:
    print("Sample size is over 10,000 - this is too large. Exiting!")
    sys.exit()
migration_matrix_filepath = str(sys.argv[4])
migration_matrix = pd.read_csv(migration_matrix_filepath, header=0, index_col=0).to_dict(orient='index')
#migration_matrix = {"PRL":{"PRL":0.999,"HMR":0.001, "LGR":0.001},
#                    "HMR":{"PRL":0.001,"HMR":0.999,"LGR":0.001},
#                    "LGR":{"PRL":0.001,"HMR":0.999,"LGR":0.001}}
#migration_matrix_filepath = "data/true_migration_prob_matrix.csv"
num_cuts = 10
ref_seq = "TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC"
# positions are 1-indexed
ref_cut_sites = [17, 43, 69, 95, 121, 147, 173, 199, 225, 251]
try:
    float(m)
    m = [float(m)] * num_cuts
except:
    if len(m.split(",")) == num_cuts:
        m = m.split(",")
        m = [float(i) for i in m]
    else:
        print("Comma-separated list of mutation rates was not of len 1 or equal to len of num_cuts. Exiting!")
        sys.exit()

## SIMULATE TREE
# patch Cassiopeia overlay function
cas.sim.Cas9LineageTracingDataSimulator.overlay_data = overlay_data_drop

# simulate a tree based on birth death model with num_extant leaves
# information on tree object: https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.data.CassiopeiaTree.html
#bd_sim = cas.sim.BirthDeathFitnessSimulator(
#    birth_waiting_distribution = lambda scale: np.random.exponential(scale),
#    initial_birth_scale = 0.5,
#    death_waiting_distribution = lambda: np.random.exponential(1.5),
#    mutation_distribution = lambda: 1 if np.random.uniform() < 0.5 else 0,
#    fitness_distribution = lambda: np.random.normal(0, .5),
#    fitness_base = 1.3,
#    num_extant = sample_num,
#    random_seed=17
#)

# run until 10,000 leaves
bd_sim = cas.sim.BirthDeathFitnessSimulator(
        #experiment_time = 280
        birth_waiting_distribution = lambda scale: np.random.exponential(scale),
        initial_birth_scale = 0.5,
        num_extant = 10000
    )
ground_truth_tree_population = bd_sim.simulate_tree()
# downsample leaves
ground_truth_tree = cas.sim.UniformLeafSubsampler(number_of_leaves=sample_num).subsample_leaves(ground_truth_tree_population)

# information on tree object: https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.data.CassiopeiaTree.html
#print("True tree:",ground_truth_tree.get_newick(record_branch_lengths=False))
#s1 = ground_truth_tree.get_newick(record_branch_lengths=True)

# overlay tissue labels for migration information
tissue_labels_df, labeled_tree = assign_tissue_labels(ground_truth_tree_population,migration_matrix)

# subsample tissue and tree migration information
tissue_labels_df, labeled_tree = subsample_tissue_labels(tissue_labels_df,ground_truth_tree)

# Write tissue labels df to output
tissue_labels_df.to_csv(outprefix + "_tissues.tsv", sep='\t', index=False)

#Write new tree with tissue labeled nodes to newick output
out_tree_tissues = outprefix + "_true_tissues.nwk"
with open(out_tree_tissues,'w') as ttt:
    ttt.write(labeled_tree.write(format=8))

# Cassipeia can take a list of values for the mutation rate, applying a different rate to each site
final_matrix = sim_chars(ground_truth_tree,m,num_cuts)
site_names = list(final_matrix.columns)
seq_arr = np.array(list(ref_seq))
# alternative way to get ref cut sites but with unequal spacing due to rounding
#ref_cut_sites = list(np.round(np.linspace(16, len(seq_arr) - 10, num_cuts)).astype(int))

## EMPIRICAL CALC EVEN CUT SITES ##
# How empirical cut sites in data were defined
# This matches the hardcoded ref_seq above in terms of RNA guide sequences
# number of base pairs buffer at start and end of amplicon
#buffer1 = 17
#buffer2 = 10
# cut_distances guaranteed to be equal to prevent biases in dropout
#cut_dist = int((len(seq_arr)-(0))/num_cuts)
# number of cut sites guaranteed to equal num_cuts
#ref_cut_sites = list(np.arange(buffer1, (len(seq_arr) - buffer2)+(buffer1+buffer2), cut_dist)[:num_cuts])

## GENERIC CALC EVEN CUT SITES ##
# Uncomment this to allow flexible number of cut sites
# This will no longer match the RNA guide sequences in ref_seq with implications for alignment
#buffer1 = 10
#buffer2 = 10
# cut_distances guaranteed to be equal to prevent biases in dropout
#cut_dist = int((len(seq_arr)-(buffer1+buffer2))/num_cuts)
# number of cut sites guaranteed to equal num_cuts
#ref_cut_sites = list(np.arange(buffer1, (len(seq_arr) - buffer2), cut_dist)[:num_cuts])

## SIMULATE MUTATIONS AND MUTANT SEQS ##
sample_names = list(final_matrix.index)
mut_dict = {}
cut_dict = {}
mut_ages = {}
row_fasta = {}
for index, row in final_matrix.iterrows():
    for site_i,s in enumerate(site_names):
        if "_" in str(row[s]):
            clean_name = int(str(row[s].split("_")[0]))
            age = float(str(row[s].split("_")[1]))
            if clean_name not in mut_ages:
                # note that age = amount of time that passed since divergence from root
                # so the oldest mutations will have the lowest numerical value
                mut_ages[clean_name] = age
            # update matrix
            final_matrix.at[index,s] = clean_name
            # update row
            row[s] = clean_name
            # might not need cut_dict but keep for now
            if s not in cut_dict:
                cut_dict[s] = ref_cut_sites[site_i]
    ages = []
    mutations = list(row)
    # get order of sites by time of mutation
    for m in mutations:
        if str(m) != '0':
            age = mut_ages[m]
        else:
            age = 0
        ages.append(age)
    # sort from oldest to youngest
    idx_by_age = np.argsort(ages)
    # keep track of deleted sites
    deleted_site_idx = []
    # mutation range
    amplicon_seq = ref_seq[:]
    amplicon_cut_sites = ref_cut_sites[:]
    for i in idx_by_age:
        if i not in deleted_site_idx and ages[i] > 0:
            mutation = mutations[i]
            if mutation not in mut_dict:
                mut_dict[mutation] = generate_indel()
            mseq = mut_dict[mutation]
            indel_len = len(mseq)
            mpos = amplicon_cut_sites[i]
            # if deletion
            if "-" in mseq:
                # add intersecting cut sites to deleted
                for i_site,s in enumerate(amplicon_cut_sites):
                    # interval comparison ignoring focal site
                    if s != mpos and mpos <= s <= (mpos+indel_len):
                        deleted_site_idx.append(i_site)
                        # update final_matrix value to missing, here encoded as -1
                        final_matrix.at[index,site_names[i_site]] = -1
                # update amplicon_seq with deletion
                amplicon_seq = mutate_seq(amplicon_seq,mpos,mseq)
            # if insertion
            else:
                amplicon_seq = mutate_seq(amplicon_seq,mpos,mseq)
                # shift amplicon cut_sites
                for i_site,s in enumerate(amplicon_cut_sites):
                    if s > mpos:
                        # increment s by indel_len
                        amplicon_cut_sites[i_site] = s + indel_len
    row_fasta[index]=amplicon_seq
# INFER TREE FROM SIMULATED MUTATIONS ## 
reconstructed_tree = cas.data.CassiopeiaTree(character_matrix = final_matrix, missing_state_indicator = -1)
greedy_solver = cas.solver.VanillaGreedySolver()
greedy_solver.solve(reconstructed_tree)

## WRITE OUTPUT FILES ##

# write trees, mutation matrix and dictionaries to output tables
out_mut = outprefix + "_mutations.tsv"
out_cut = outprefix + "_cut_positions.tsv"
out_matrix = outprefix + "_indel_character_matrix.tsv"
with open(out_mut,'a') as outm:
    for k,v in mut_dict.items():
        age = mut_ages[k]
        outline = str(k) + "\t" + str(v) + "\t" + str(age) + "\n"
        outm.write(outline)
with open(out_cut,'a') as outc:
    for k,v in cut_dict.items():
        outline = str(k) + "\t" + str(v) + "\n"
        outc.write(outline)
final_matrix.to_csv(out_matrix, sep='\t')

# inferred tree
out_tree_true = outprefix + "_true.nwk"
out_tree_infer = outprefix + "_inferred.nwk"
with open(out_tree_true,'w') as tt:
    tt.write(ground_truth_tree.get_newick(record_branch_lengths=False))
with open(out_tree_infer,'w') as it:
    it.write(reconstructed_tree.get_newick())

# Fasta sequence for each simulated sample
outfasta = outprefix + ".fa"
with open(outfasta,'a') as o:
    outstr = ">ref" + "\n" + ref_seq + "\n"
    o.write(outstr)

# Write tissue specific fasta files
tissues = tissue_labels_df['tissue'].unique()
for tis in tissues:
    globals()[f"outfasta_{tis}"] = outprefix + "_" + tis + ".fa"
    with open(globals()[f"outfasta_{tis}"],'a') as o:
        outstr = ">ref" + "\n" + ref_seq + "\n"
        o.write(outstr)

for i,seq in row_fasta.items():
    seq = seq.replace('-', '')
    outstr = ">" + str(i) + "\n" + seq + "\n"
    with open(outfasta,'a') as o:
        o.write(outstr)
    
    tissue = tissue_labels_df.loc[tissue_labels_df['node'] == i, 'tissue'].values[0]
    with open(globals()[f"outfasta_{tissue}"],'a') as t:
        t.write(outstr)
