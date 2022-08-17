import sys
import math
import pandas as pd
import ete3
from ete3 import Tree
from Bio import Phylo
import numpy as np
from scipy.optimize import minimize

#read in Phylogenetic Tree
#tree = Phylo.read(sys.argv[1], 'newick')
asv_tab = sys.argv[2]
cp_group = sys.argv[3]

out_tree = sys.argv[3] + "_tree.txt"
out_labels = sys.argv[3] + "_labels.txt"
out_colors = sys.argv[3] + "_colors.txt"

# read in nested ASV dict
#asv_names,sample,group
#ASV001,LGR,CP01
asv_dict = {}
# group is first key, asv is second key
with open(asv_tab,'r') as tab:
    for l in tab:
        l = l.strip()
        l = l.split(",")
        if l[2] not in asv_dict:
            asv_dict[l[2]] = {}
            asv_dict[l[2]][l[0]] = [l[1]]
        else:
            if l[0] in asv_dict[l[2]]:
                if l[1] not in asv_dict[l[2]][l[0]]:
                    asv_dict[l[2]][l[0]].append(l[1])
            else:
                asv_dict[l[2]][l[0]] = [l[1]]

# leaves in CP group
cp_dict = asv_dict[cp_group]
print(cp_dict)
labelfile = open(out_labels, 'w')
tissues = []
for k,v in cp_dict.items():
    for tissue in v:
        if tissue not in tissues:
            tissues.append(tissue)
        outline = str(k) + " " + str(tissue)
        labelfile.write("%s\n" % outline)
labelfile.close()

colorfile = open(out_colors, 'w')
for i,t in enumerate(tissues):
    outline = str(t) + " " + str(i+1)
    colorfile.write("%s\n" % outline)
colorfile.close()

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
    
