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
            print(node,tree.get_time(node),character_array, cuts_remaining)
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





sample_num = 10
bd_sim = cas.sim.BirthDeathFitnessSimulator(
    birth_waiting_distribution = lambda scale: np.random.exponential(scale),
    initial_birth_scale = 0.5,
    death_waiting_distribution = lambda: np.random.exponential(1.5),
    mutation_distribution = lambda: 1 if np.random.uniform() < 0.5 else 0,
    fitness_distribution = lambda: np.random.normal(0, .5),
    fitness_base = 1.3,
    num_extant = sample_num,
    random_seed=17
)
ground_truth_tree = bd_sim.simulate_tree()
np.random.seed(seed=None)
mut_rate = 0.1
cas.sim.Cas9LineageTracingDataSimulator.overlay_data = overlay_data_drop

#mytest1 = cas.sim.Cas9LineageTracingDataSimulator()
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
#lt_sim.overlay_data(ground_truth_tree)
#mytest1 = Cas9LineageTracingDataSimulator()
lt_sim.overlay_data(ground_truth_tree)

df = ground_truth_tree.character_matrix
print(df)
site_names = list(df.columns)

for index, row in df.iterrows():
    #print(index,list(row))
    print(list(row))
    for s in site_names:
        print(index,s,row[s])
        if "_" in str(row[s]):
            clean_name = int(str(row[s].split("_")[0]))
            df.at[index,s] = clean_name
print(df)
print(ground_truth_tree.get_newick())

