import sys
from ete3 import Tree

def get_tree_depth(tree):
    max_depth = 0
    for l in [leaf.name for leaf in tree.iter_leaves()]:
        leaf_depth = -1
        node = tree.search_nodes(name=l)[0]
        while node:
            node = node.up
            leaf_depth += 1
        if leaf_depth > max_depth:
            max_depth = leaf_depth
    return max_depth

primary_tissue = sys.argv[3]

labels = {}
with open(sys.argv[1],'r') as lab:
    for l in lab:
        l = l.strip().split(' ')
        node = l[0]
        tissue =  l[1]
        if node not in labels:
            labels[node] = tissue
tabular_tree = []
with open(sys.argv[2],'r') as table:
    for l in table:
        l = l.strip().split(' ')
        parent = l[0]
        child = l[1]
        distance = 1
        tabular_tree.append((parent,child,distance))

tree = Tree.from_parent_child_table(tabular_tree)
tree_depth = get_tree_depth(tree)
print(tree)

if tree_depth > 1:
    # ignore root node
    for node in tree.iter_descendants("levelorder"):
        if not node.is_leaf():
            # get all child nodes
            node_descendants = []
            node_node_descendants = []
            for nd in node.children:
                if not nd.is_leaf():
                    node_node_descendants.append(nd.name)
                else:
                    node_descendants.append(nd.name)
            if len(node_node_descendants) > 0:
                print(node.name,node)

            # no child nodes belong to primary tissue
            #print(node)
            # 
        #node_tissue = labels[node.name]
        #print(node.name,node_tissue)

