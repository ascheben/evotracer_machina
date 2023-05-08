#!/usr/bin/env python3
import sys
import csv
import numpy as np
import cassiopeia as cas
from cassiopeia.data import CassiopeiaTree

def all_children(tree,node):
    if not tree.children(node):
        return []
    result = [ (list(tree.children(node))) ]
    for child in tree.children(node):
        result.extend(all_children(tree,child))
    return result

def num_mutated_leaves(cas_tree,mut_prob):
    tree = cas_tree
    mutated_leaves = set()
    tis_labels = {}
    root = tree.root
    # Traverse the tree and allow each node to mutate
    for node in tree.nodes:
        if node not in tis_labels:
            tis_labels[node] = "wt"
        if node == root:
            continue
        if tree.is_leaf(node) == False and tis_labels[node] == "wt":
            node_children = all_children(tree,node)
            node_children = [item for sublist in node_children for item in sublist]
            node_children_tip = []
            for nc in node_children:
                if tree.is_leaf(nc):
                    node_children_tip.append(nc)
            genotype = np.random.choice(["mut","wt"], p=[mut_prob,1-mut_prob])
            tis_labels[node] = genotype
            if tis_labels[node] == "mut":
                for c in node_children:
                    if tree.is_leaf(c):
                        mutated_leaves.add(c)
                    else:
                        tis_labels[node] = "mut"

    return len(mutated_leaves)
## PARAMETERS ##
mut_prob = float(sys.argv[1])
output_name=sys.argv[2]
sample_num = 100
try:
    bd_sim = cas.sim.BirthDeathFitnessSimulator(
        #experiment_time = 280
        birth_waiting_distribution = lambda scale: np.random.exponential(scale),
        initial_birth_scale = 0.5,
        num_extant = 10000
    )
    ground_truth_tree = bd_sim.simulate_tree()
    # downsample
    ground_truth_tree = cas.sim.UniformLeafSubsampler(number_of_leaves=sample_num).subsample_leaves(ground_truth_tree)
    #print(ground_truth_tree.get_newick(record_branch_lengths=True))
    mut = num_mutated_leaves(ground_truth_tree,mut_prob)
    freq_mut = mut / sample_num
    print(mut_prob,mut,freq_mut)
    #write output to a csv file to simplify running in parallel
    data = [output_name ,mut_prob, mut, sample_num, freq_mut]
    filename=f"{output_name}_freq_mut_from_rate_data.csv"
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data)
except:
    print("Simulation failed. Try again!")
