import sys
import ete3
from ete3 import Tree

t = Tree(sys.argv[1])
all_leaves = []
for leaf in t:
    all_leaves.append(leaf.name)

leaf1 = "n83"
leaf2 = "n27"

A = t&leaf1
B = t&leaf2
tmrca = A.get_distance(B) / 2
#read in Phylogenetic Tree
print("TMRCA age for",leaf1," and ", leaf2, " is ", tmrca)

