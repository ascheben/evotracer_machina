import sys
import re
from ete3 import Tree

#grep "\->" T-PRL-S-binarized.dot| cut -d' ' -f1,3|sed 's/^\t//' > T-PRL-S-binarized_parent_child.txt
#grep label T-PRL-S-binarized.dot| tr ',' ' '| cut -d' ' -f1,4,5| sed 's/^\t//'| sed 's/color=//'| sed 's/label="//'| sed 's/"//g'| sed 's/]//'| sed 's/\\n/__/' |sed 's/^\t//'> T-PRL-S-binarized_nodes.txt

nodes = sys.argv[1]
edges = sys.argv[2]
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
        # ignore redundant labels of false polytomy resolutions
        #if name.count('^') < 2 or name.split("^")[2]=="0":
            #print("label",col,name)
            
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

        # ignore redundant parent-child relationships caused by false polytomy resolutions
        #if not parent_name.count('^') > 1 and not child_name.count('^') > 1:
        #    print("tree",parent_name,child_name)
        #elif parent_name.count('^') > 1 and not child_name.count('^') > 1:
        #    if not parent_name.split("^")[2]=="0":
        #        parent_name = parent_name.split("^")[0] + "^" + parent_name.split("^")[1] + "^0"
        #    print("tree",parent_name,child_name)
        #elif not parent_name.count('^') > 1 and child_name.count('^') > 1:
        #    if child_name.split("^")[2]=="0":
        #        print("tree",parent_name,child_name)


# traverse parent node
# simplify graph
# if a parent A has tissue_label X, and it's parent B has label X, and it's children both have X
# then delete the parent A and assign the children to B

# Add correct node labels to tree
#for i,node in enumerate(tree.traverse("preorder")):
#    if node.is_leaf():
#        ntype = "leaf"
#    else:
#        ntype = "internal"




# init non-redundant parent dict
nr_pardict = {}
redundant_nodes = []
nr_nodes = []
for p,c in pardict.items():
    if p not in chidict:
        # is root node
        parent_parent = None
    else:
        parent_label = coldict[nodecoldict[p]]
        parent_parent = chidict[p]
        parent_parent_label = coldict[nodecoldict[parent_parent]]
        children_label = set()
        for child in c:
            children_label.add(coldict[nodecoldict[child]])
        children_label = list(children_label)
        if len(children_label) == 1:
            if children_label[0] == parent_parent_label and children_label[0] == parent_label: 
                #print("Redundant node - Parent node, parent_parent, children, parent_parent_label,children_label:", p, parent_parent, c, parent_parent_label,children_label[0])
                # remove redundant
                redundant_nodes.append(p)
                if parent_parent in nr_pardict:
                    cur_children = nr_pardict[parent_parent]
                    # remove parent
                    cur_children = [s for s in cur_children if s != p]
                    # add children
                    cur_children = cur_children + c
                    nr_pardict[parent_parent] = cur_children
                else:
                    nr_pardict[parent_parent] = c
            else:
                #print("Non-Redundant node - Parent node, parent_parent, children, parent_parent_label,children_label:", p, parent_parent, c, parent_parent_label,children_label[0])
                nr_pardict[p] = c
                nr_nodes.append(p)


tree = Tree.from_parent_child_table(tabular_tree)
for r in redundant_nodes:
    badnode = tree.search_nodes(name=r)[0]
    badnode.delete()

for node in tree.traverse("preorder"):
    if node.is_leaf():
        if node.name not in nr_nodes:
            nr_nodes.append(node.name)

#flatlist=[element for sublist in list(nr_pardict.values()) for element in sublist]
#nodes = list(set(list(nr_pardict.keys()) + flatlist))
for n in nr_nodes:
    ncol = nodecoldict[n]
    print("label",ncol,n)

for node in tree.traverse("preorder"):
    if not node.is_leaf():
        parent_name = node.name
        for c in node.children:
            child_name = c.name
            print("tree",parent_name,child_name)

#print(tree)
#for p,c in nr_pardict.items():
#    for child in c:
#        print("tree",p,child)
# nodes
#277 3 ASV030_PRL\nPRL
#276 1 ASV222_LVM\nLVM

# parent child
#274 276
#274 275
