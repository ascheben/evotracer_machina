#!/usr/bin/env python3
import sys
import os
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

out_tree = sys.argv[3] + "_tree_split.txt"
out_labels = sys.argv[3] + "_labels_split.txt"
out_colors = sys.argv[3] + "_colors.txt"
out_newick = sys.argv[3] + ".newick"
out_failed = "FailedCP.txt"

# hardcode primary tissue to add dummy where its missing
#primary_tissue = "PRL"
primary_tissue = sys.argv[4]

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
#print(cp_dict)
primary_absent = True
labelfile = open(out_labels, 'w')
tissues = []
split_samples = {}
for k,v in cp_dict.items():
    tissue_count = len(v)
    for tissue in v:
        if tissue == primary_tissue:
            primary_absent = False
        if tissue not in tissues:
            tissues.append(tissue)
        if tissue_count == 1:
            outline = str(k) + " " + str(tissue)
        else:
            outline = str(k) + "_" + str(tissue) + " " + str(tissue)
            if k not in split_samples:
                split_samples[k] = [str(k) + "_" + str(tissue)]
            else:
                split_samples[k].append(str(k) + "_" + str(tissue))
        labelfile.write("%s\n" % outline)
if primary_absent:
    outline = "ASVXXX" + " " + primary_tissue
    labelfile.write("%s\n" % outline)
labelfile.close()

if primary_absent:
    tissues.append(primary_tissue)

colorfile = open(out_colors, 'w')
for i,t in enumerate(tissues):
    outline = str(t) + " " + str(i+1)
    colorfile.write("%s\n" % outline)
colorfile.close()

species_list = list(cp_dict.keys())
prune_tree = Tree(sys.argv[1]) 
prune_tree.prune(species_list,preserve_branch_length=True)
prune_tree.write(format=1, outfile=out_newick)

#tree = Phylo.read(sys.argv[1], 'newick')
tree = prune_tree

#common_ancestor = tree.common_ancestor(species_list)
# if tree tips dont match expected due to unexpected topology delete all output and exit
#observed_species_list = []
#for t in common_ancestor.get_terminals():
#    observed_species_list.append(t.name)
#if set(species_list) == set(observed_species_list):
#    Phylo.write(common_ancestor, out_newick, "newick")
#    tree = Tree(out_newick)
#else:
#    overlap = len(set(species_list).intersection(set(observed_species_list)))
#    print(cp_group,"original ASV names was ", len(species_list), "new ASV names was ", len(observed_species_list),"overlap was ", overlap)
    #print("No match between",species_list," and ",observed_species_list)
#    os.remove(out_labels)
#    os.remove(out_colors)
#    cp_group
#    with open(out_failed, "a") as myfile:
#            myfile.write(cp_group + "\n")
#    sys.exit()


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
        if c.name in split_samples:
            for tissue in split_samples[c.name]:
                outline = str(n.name) + " " + str(tissue)
                treefile.write("%s\n" % outline)
        else:
            outline = str(n.name) + " " + str(c.name)
            treefile.write("%s\n" % outline)
if primary_absent:
    outline = "0" + " " + "ASVXXX"
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
    
