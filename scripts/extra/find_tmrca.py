import sys
import math
import pandas as pd
import ete3
from ete3 import Tree
from Bio import Phylo
import numpy as np
from scipy.optimize import minimize

species_list = list(cp_dict.keys())
tree = Phylo.read(sys.argv[1], 'newick')
common_ancestor = tree.common_ancestor(species_list)
Phylo.write(common_ancestor, "subtree.newick", "newick")
tree = Tree("subtree.newick")

#tree = Tree(sys.argv[1])

node2leaves = tree.get_cached_content()
nodes = {}
node_name_counter = 0
leaves = []
for leaf in tree.iter_leaves():
    leaves.append(leaf.name)
    # read feature from file
    #leaf_features = asv_dict[leaf.name] 

for n in tree.traverse():
    if n.name == "":
        n.name = node_name_counter
        node_name_counter += 1

treefile = open(out_tree, 'w')

for n in tree.traverse("postorder"):
    #print("node %s contains %s tips" %(n.name, len(node2leaves[n])))
    #print("Node children:")
    # assign most common features in children
    # for each feature, collect the most frequent in the children
    #print(n.name)
    for c in n.children:
        outline = str(n.name) + " " + str(c.name)
        treefile.write("%s\n" % outline)
treefile.close()
    #if len(n.children) == 0:
    #    pass
        # already in asv_dict
        #node_features = asv_dict[n.name]
    #else:
    #    asv_dict[n.name] = node_features
    #    for feat in node_features.keys():
    #        child_features = []
    #        for child in n.children:
    #            child_features.append(asv_dict[child.name][feat])
    #        if -1 in child_features:
    #            common = -1
    #        else:
    #            common = most_common(child_features)
    #        asv_dict[n.name][feat] = common
    
