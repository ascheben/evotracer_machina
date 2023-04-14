#!/usr/bin/env python3
import sys
import cProfile
from collections import defaultdict
import copy
import numpy as np
import pandas as pd
import time
import random
import cassiopeia as cas
from ete3 import Tree

#import pdb
#import networkx as nx
#from tqdm.auto import tqdm
#import matplotlib.pyplot as plt
#from cassiopeia.solver import missing_data_methods
#import dendropy
#from dendropy.calculate import treecompare
#from dendropy.simulate import treesim

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

def generate_indel(max_size):
    '''
    Take sequence and position and randomly introduce an
    indel at the given position. Returns the mutated sequence.
    '''
    del_prob = 0.6
    ins_prob = 0.4
    #max_size = 15
    #max_size = 3
    mut_size = int(np.random.exponential(6,1)[0])
    while mut_size == 0 or mut_size > max_size:
        mut_size = int(np.random.exponential(6,1)[0])
    mut_type = random.choices(population=['ins','mut'],weights=[ins_prob,del_prob],k=1)[0]
    #mut_size = random.randint(1, 10)
    if mut_type == 'ins':
        #generate random kmer of mut_size bp
        mut_bases = ''.join(random.choices(['A','C','T','G'], k=mut_size))
    else:
        mut_bases = mut_size * "-"
    return mut_bases

def sim_chars(tree,mut_rate):
    np.random.seed(seed=None)
    lt_sim = cas.sim.Cas9LineageTracingDataSimulator(
        number_of_cassettes = 5,
        size_of_cassette = 1,
        mutation_rate = mut_rate,
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
            tis_labels[node] = np.random.choice(tissue_array, p=prob_vec)
        # Label leaves not by probability but by their parent
        if tree.is_leaf(node) == True:
            tis_labels[node] = prev_tissue
    
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

    return tissues_df, tree_labeled
    
    
outprefix = sys.argv[1]
high_mut_rate = float(sys.argv[2])
low_mut_rate = float(sys.argv[3])
max_size = int(sys.argv[4])
sample_num = int(sys.argv[5])
migration_matrix_filepath = str(sys.argv[6])

#max_size = 3
#migration_matrix = {"prostate":{"prostate":0.34,"lung":0.33, "liver":0.33},
#                    "lung":{"prostate":0.33,"lung":0.34,"liver":0.33},
#                    "liver":{"prostate":0.33,"lung":0.33,"liver":0.34}}
#migration_matrix_filepath = "data/true_migration_prob_matrix.csv"

migration_matrix = pd.read_csv(migration_matrix_filepath, header=0, index_col=0).to_dict(orient='index')

# simulate a tree based on birth death model with num_extant leaves
bd_sim = cas.sim.BirthDeathFitnessSimulator(
    birth_waiting_distribution = lambda scale: np.random.exponential(scale),
    initial_birth_scale = 0.5,
    death_waiting_distribution = lambda: np.random.exponential(1.5),
    mutation_distribution = lambda: 1 if np.random.uniform() < 0.5 else 0,
    fitness_distribution = lambda: np.random.normal(0, .5),
    fitness_base = 1.3,
    num_extant = sample_num,
    #num_extant = 10,
    random_seed=17
)
ground_truth_tree = bd_sim.simulate_tree()
# information on tree object: https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.data.CassiopeiaTree.html
#print("True tree:",ground_truth_tree.get_newick(record_branch_lengths=False))
#s1 = ground_truth_tree.get_newick(record_branch_lengths=True)

# overlay tissue labels for migration information
tissue_labels_df, labeled_tree = assign_tissue_labels(ground_truth_tree,migration_matrix)

# Write tissue labels df to output
tissue_labels_df.to_csv(outprefix + "_tissues.tsv", sep='\t', index=False)

#Write new tree with tissue labeled nodes to newick output
out_tree_tissues = outprefix + "_true_tissues.nwk"
with open(out_tree_tissues,'w') as ttt:
    ttt.write(labeled_tree.get_newick())



#mut_rates = [0.1,0.01]
mut_rates = [high_mut_rate,low_mut_rate]
for i,m in enumerate(mut_rates):
    if i != 0:
        suffix = "_" + str(i) + "_" + str(m)
        character_matrix = character_matrix.join(sim_chars(ground_truth_tree,m),rsuffix=suffix)
    else:
        character_matrix = sim_chars(ground_truth_tree,m) 
# join the matrices on the index col that is the tree leaf label
final_matrix = character_matrix
reconstructed_tree = cas.data.CassiopeiaTree(character_matrix = final_matrix, missing_state_indicator = -1)
greedy_solver = cas.solver.VanillaGreedySolver()
greedy_solver.solve(reconstructed_tree)

# Get the reconstructed tree newick
#s2 = reconstructed_tree.get_newick() 

# Resolve polytomy for RF calculation
#t = Tree(reconstructed_tree.get_newick())
#t.resolve_polytomy(recursive=True)
#s2 = t.write()
#nj_solver = cas.solver.NeighborJoiningSolver(add_root=False)

# Calculate RF distance between true and reconstructed tree
# https://dendropy.org/primer/treestats.html
# establish common taxon namespace
#tns = dendropy.TaxonNamespace()
#taxa = dendropy.TaxonNamespace(reconstructed_tree.leaves)

# ensure all trees loaded use common namespace
#tree1 = dendropy.Tree.get(
#        data=s1,
#        schema='newick',
#        taxon_namespace=taxa)
#tree2 = dendropy.Tree.get(
#        data=s2,
#        schema='newick',
#        taxon_namespace=taxa)
#ntax = len(reconstructed_tree.leaves)
#t = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=ntax)


# prepare mutated sequence
ref_seq = "TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC"
# positions are 1-indexed
ref_cut_sites = [17, 43, 69, 95, 121, 147, 173, 199, 225, 251] 
#ref_border_sites = [1, 26, 52, 78, 104, 130, 156, 182, 208, 234]
site_names = list(final_matrix.columns)
sample_names = list(final_matrix.index)
site_indels = set()
mut_dict = {}
cut_dict = {}

for m,s in enumerate(site_names):
    cut_dict[s] = ref_cut_sites[m]
    indel_list = final_matrix[s].tolist()
    mutations = []
    #print("checking site:",m,s)
    for i in indel_list:
        #print("checking indel:", i)
        if str(i) != '0':
            site_indels.add(i)
            #print("Added to site_indel set:", site_indels)
            if i not in mut_dict:
                indel_seq = generate_indel(max_size)
                #print("Indel not in mut_dict - generated indel:",indel_seq)
                # ensure all indels at each site are unique
                while indel_seq in mutations:
                    #print("Indel seq",indel_seq,"already exists at site:",mutations)
                    indel_seq = generate_indel(max_size)
                #print("Adding new indel to site",indel_seq)
                mutations.append(indel_seq)
                mut_dict[i] = indel_seq

# write trees, mutation matrix and dictionaries to output tables
out_mut = outprefix + "_mutations.tsv"
out_cut = outprefix + "_cut_positions.tsv"
out_matrix = outprefix + "_indel_character_matrix.tsv"
with open(out_mut,'a') as outm:
    for k,v in mut_dict.items():
        outline = str(k) + "\t" + str(v) + "\n"
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

# Write overall fasta file
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


for index, row in final_matrix.iterrows():
    row_fasta = ref_seq
    tissue = tissue_labels_df.loc[tissue_labels_df['node'] == index, 'tissue'].values[0]
    #print("refseq",ref_seq)
    #print(index,row)
    insertion_pos_dict = {}
    #print("ins_pos_dict",insertion_pos_dict)
    for s in site_names:
        indel_id = row[s]
        cut_site = cut_dict[s]
        # if not no indel
        if indel_id != 0:
            #print("Mutation at site")
            mut_bases = mut_dict[indel_id]
            #print("mut_bases:",mut_bases)
            insertion_shift = 0
            for k,v in insertion_pos_dict.items():
                # preceding insertion
                if k < cut_site:
                    insertion_shift += len(v)
            #print("Final insertion shift:",insertion_shift)
            cut_site += insertion_shift
            row_fasta = mutate_seq(row_fasta,cut_site,mut_bases)
            #print("cut site,rowfasta", cut_site,row_fasta)
            # track insertions
            if '-' not in mut_bases:
                insertion_pos_dict[cut_site] = mut_bases
                #print("Added cut site to ins_pos_dict", cut_site, mut_bases)
    row_fasta = row_fasta.replace('-', '')
    outstr = ">" + index + "\n" + row_fasta + "\n"
    with open(outfasta,'a') as o:
        o.write(outstr)
        
    with open(globals()[f"outfasta_{tissue}"],'a') as t:
        t.write(outstr)
