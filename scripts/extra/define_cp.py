import sys
from ete3 import Tree

def find_cp_in_tree(tree):
    #print(tabular_tree)
    #print(tree.write())
    # ignore root node
    for node in tree.iter_descendants("levelorder"):
        if not node.is_leaf() and node.up.is_root():
#            sisters = []
#            sister_descendants = []
#            sister_descendant_freqs = []
#            for n in node.get_sisters():
#                sisters.append(n.name)
#                descendants = n.get_leaf_names()
#                sister_descendants = sister_descendants + descendants
#                for d in descendants:
#                    if d in freq_dict[cp]:
#                        sister_descendant_freqs.append(freq_dict[cp][d])
#            node_descendants = []
#            node_descendant_freqs = []
#            for nd in node.get_leaf_names():
#                node_descendants.append(nd)
#                if nd in freq_dict[cp]:
 


t = Tree(sys.argv[1])

find_cp_in_tree(t)
