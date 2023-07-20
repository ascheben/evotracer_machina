import sys
import re
from ete3 import Tree

#grep "\->" T-PRL-S-binarized.dot| cut -d' ' -f1,3|sed 's/^\t//' > T-PRL-S-binarized_parent_child.txt
#grep label T-PRL-S-binarized.dot| tr ',' ' '| cut -d' ' -f1,4,5| sed 's/^\t//'| sed 's/color=//'| sed 's/label="//'| sed 's/"//g'| sed 's/]//'| sed 's/\\n/__/' |sed 's/^\t//'> T-PRL-S-binarized_nodes.txt

# Revise this method:
# 1) mark as safe all nodes that exist in the original tree based on exact set of children
# 2) Remove spurious self-transitions
#- If parent of a node is same:
#  delete that node unless:
#      a) its constrained by mutation-supported lineage in CP03_tree_split.txt
#      b) its a migration node (its parent is diff color)

nodes = sys.argv[1]
edges = sys.argv[2]
split_tree = sys.argv[3]
coldict = {}
nodedict = {}
nodecoldict = {}

with open(nodes,'r') as n:
    for l in n:
        l = l.strip()
        l = l.split(" ")
        name = l[2]
        col = l[1]
        nodeid = l[0]
        if "ASV" in name:
            split_name = re.split('__',name)
            tissue = split_name[1]
            name = split_name[0]
            if col not in coldict:
                coldict[col] = tissue
        else:
            # assumes coldict has this value
            try:
                tissue = coldict[col]
            except:
                tissue = "NA"

        if name not in nodecoldict:
            nodecoldict[name] = col
        if nodeid not in nodedict:
            nodedict[nodeid] = name
            
pardict = {}
chidict = {}
tabular_tree = []

with open(edges,'r') as e:
    for l in e:
        l = l.strip()
        l = l.split(" ")
        parent = l[0]
        child = l[1]
        parent_name = nodedict[parent]
        child_name = nodedict[child]
        tabular_tree.append((parent_name,child_name,1))
        if parent_name not in pardict:
            pardict[parent_name] = [child_name]
        else:
            pardict[parent_name].append(child_name)
        if child_name not in chidict:
            chidict[child_name] = parent_name

s_tabular_tree = []
with open(split_tree,'r') as s:
    for l in s:
        l = l.strip().split(" ")
        parent = l[0]
        child = l[1]
        s_tabular_tree.append((parent,child,1))

redundant_nodes = []
tree = Tree.from_parent_child_table(tabular_tree)
# split cassiopeia tree
stree = Tree.from_parent_child_table(s_tabular_tree)
conf_nodes = []
for node in stree.iter_descendants("preorder"):
    if not node.is_leaf():
        descendants = []
        # get list of descendants (tips only)
        for d in node.iter_descendants("preorder"):
            if d.is_leaf():
                descendants.append(d.name)
        conf_nodes.append(descendants)
conf_nodes_found = []
# skip root
for node in tree.iter_descendants("preorder"):
    if not node.is_leaf():
        descendants = []
        # get list of descendants (tips only)
        for d in node.iter_descendants("preorder"):
            if d.is_leaf():
                descendants.append(d.name)
        # is this node supported by the raw MP tree
        node_is_supported = False
        for c in conf_nodes:
            if sorted(c) == sorted(descendants):
                node_is_supported = True
        # only allow one node per supported lineage
        for f in conf_nodes_found:
            if sorted(f) == sorted(descendants):
                node_is_supported = False
        if node_is_supported:
            conf_nodes_found.append(descendants)
        else:
            # check if spurious self transition node
            parent_node = node.up.name
            parent_node_tissue = coldict[nodecoldict[parent_node]]
            node_tissue = coldict[nodecoldict[node.name]]
            if node_tissue == parent_node_tissue:
                redundant_nodes.append(node.name)

for r in redundant_nodes:
    badnode = tree.search_nodes(name=r)[0]
    badnode.delete()

for node in tree.traverse("preorder"):
    ncol = nodecoldict[node.name]
    print("label",ncol,node.name)

for node in tree.traverse("preorder"):
    if not node.is_leaf():
        parent_name = node.name
        for c in node.children:
            child_name = c.name
            print("tree",parent_name,child_name)

