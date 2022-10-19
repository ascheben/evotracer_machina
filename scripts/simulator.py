import cProfile
from collections import defaultdict
import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import time
from tqdm.auto import tqdm

import cassiopeia as cas
from cassiopeia.solver import missing_data_methods
import dendropy
from dendropy.calculate import treecompare
from dendropy.simulate import treesim
from ete3 import Tree

def sim_chars(tree,mut_rate):
    np.random.seed(seed=None)
    lt_sim = cas.sim.Cas9LineageTracingDataSimulator(
        number_of_cassettes = 5,
        size_of_cassette = 1,
        mutation_rate = mut_rate,
        state_generating_distribution = lambda: np.random.exponential(1e-5),
        number_of_states = 10,
        state_priors = None,
        heritable_silencing_rate = 9e-4,
        stochastic_silencing_rate = 0.1,
        heritable_missing_data_state = -1,
        stochastic_missing_data_state = -1,
    )
    lt_sim.overlay_data(tree)
    character_matrix = ground_truth_tree.character_matrix
    return character_matrix
    



# simulate a tree based on birth death model with num_extant leaves
bd_sim = cas.sim.BirthDeathFitnessSimulator(
    birth_waiting_distribution = lambda scale: np.random.exponential(scale),
    initial_birth_scale = 0.5,
    death_waiting_distribution = lambda: np.random.exponential(1.5),
    mutation_distribution = lambda: 1 if np.random.uniform() < 0.5 else 0,
    fitness_distribution = lambda: np.random.normal(0, .5),
    fitness_base = 1.3,
    num_extant = 20,
    random_seed=17
)
ground_truth_tree = bd_sim.simulate_tree()
# information on tree object: https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.data.CassiopeiaTree.html
print("True tree:",ground_truth_tree.get_newick(record_branch_lengths=False))
s1 = ground_truth_tree.get_newick(record_branch_lengths=True)

mut_rates = [0.1,0.5]
for i,m in enumerate(mut_rates):
    if i != 0:
        suffix = "_" + str(i) + "_" + str(m)
        character_matrix = character_matrix.join(sim_chars(ground_truth_tree,m),rsuffix=suffix)
    else:
        character_matrix = sim_chars(ground_truth_tree,m) 
# join the matrices on the index col that is the tree leaf label
final_matrix = character_matrix
#final_matrix = ground_truth_tree.character_matrix.join(ground_truth_tree.character_matrix,rsuffix='_right')
# rsuffix is necessary to ensure unique column names

# the character_matrix is a pandas dataframe
# this means that we can simulate multiple dataframes with variable cutting rates
# and then join them together


reconstructed_tree = cas.data.CassiopeiaTree(character_matrix = final_matrix, missing_state_indicator = -1)

greedy_solver = cas.solver.VanillaGreedySolver()
greedy_solver.solve(reconstructed_tree)
# Get the reconstructed tree newick
print("Reconstructed:",reconstructed_tree.get_newick())
#s2 = reconstructed_tree.get_newick() 

# Resolve polytomy for RF calculation
t = Tree(reconstructed_tree.get_newick())
#t = Tree("(( (a, b, c), (d, e, f, g)), (f, i, h));")
t.resolve_polytomy(recursive=True)
#print("Clean:", t)
s2 = t.write()
print("Bifurcating reconstructed:",s2)
# make NJ tree
nj_solver = cas.solver.NeighborJoiningSolver(add_root=False)
# NJ solver wont work without explicit root
#nj_solver.solve(reconstructed_tree)

# Calculate RF distance between true and reconstructed tree
# https://dendropy.org/primer/treestats.html
#s1 = "(a,(b,(c,d)));"
#s2 = "(a,(d,(b,c)));"

# establish common taxon namespace
tns = dendropy.TaxonNamespace()
taxa = dendropy.TaxonNamespace(reconstructed_tree.leaves)

# ensure all trees loaded use common namespace
tree1 = dendropy.Tree.get(
        data=s1,
        schema='newick',
        taxon_namespace=taxa)
tree2 = dendropy.Tree.get(
        data=s2,
        schema='newick',
        taxon_namespace=taxa)

#print(reconstructed_tree.leaves)
#taxa = dendropy.TaxonNamespace(reconstructed_tree.leaves)
ntax = len(reconstructed_tree.leaves)
t = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=ntax)
print(t.as_string(schema='newick',suppress_edge_lengths=True))
## Unweighted Robinson-Foulds distance
print(treecompare.symmetric_difference(tree1, tree2))
print(treecompare.symmetric_difference(tree1, tree1))
print("Null:",treecompare.symmetric_difference(tree1, t))

