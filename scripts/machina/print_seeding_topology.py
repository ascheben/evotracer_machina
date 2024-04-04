import sys
from collections import Counter
from ete3 import Tree
from itertools import groupby

def get_index(mylist,myname):
    for i,n in enumerate(mylist):
        if myname in n:
            return i

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

def find_ABA(seeding_path):
    '''
    find patterns of tissue transitions indicating
    bidirectional seeding or reseeding
    '''
    ABA = []
    # collapse redundancy in seeding path
    seeding_path = [key for key, _group in groupby(seeding_path)]
    for i,s in enumerate(seeding_path):
        try:
            second = seeding_path[i+1]
            third = seeding_path[i+2]
            if s == third and s != second:
                pattern = [s,second,third]
                ABA.append(pattern)
        except:
            pass
    return ABA

def seeding_topology_tree(tabular_tree,tissue_dict,ptissue):
    tree = Tree.from_parent_child_table(tabular_tree)
    #tree_depth = get_tree_depth(tree)
    parallel_seeding_edges = {}
    seeding_type = []
    edges = []
    #print(tree)
    #print(tissue_dict)
    # root is skipped but we assume it is from prostate
    #primary = "PRL"
    primary = ptissue
    for i,node in enumerate(tree.traverse("preorder")):
        if not node.is_leaf():
            #print(i,node.name,node.children)
            parent_tis = tissue_dict[node.name]
            children_tis = []
            children_list = []
            primary_children = []
            primary_children_tis = []
            for n,child in enumerate(node.children):
                child_tis = tissue_dict[child.name]
                if child_tis != parent_tis:
                    children_tis.append(child_tis)
                    children_list.append(child.name)
                else:
                    primary_children.append(child.name)
                    primary_children_tis.append(child_tis)
                #print(i,node.name,child.name,parent_tis,child_tis)
            # is there parallel seeding?
            if len(set(children_tis))>1:
                # find the edges to exclude
                # here all non-primary children seeded in parallel are counted as one event
                parallel_seeding_edges[node.name] = children_list
                if parent_tis == primary:
                    seeding_type.append("primary_parallel_seeding")
                else:
                    seeding_type.append("metastatic_parallel_seeding")
                for c in primary_children:
                    edges.append([parent_tis,parent_tis])
            else:
                for c in primary_children_tis:
                    edges.append([parent_tis,c])
                for c in children_tis:
                    edges.append([parent_tis,c])
    # asv tissue for primary and cascade
    #print("New edges",cur_cp,edges)
    for pair in edges:
        if len(pair) == 2:
            if pair[0] == primary and pair[1] != primary:
                #pri2met_tissues.add(pair[1])
                seeding_type.append("primary_mono_seeding")
            elif pair[0] == primary and pair[1] == primary:
                seeding_type.append("primary_confined")
            elif pair[0] != primary and pair[1] != primary and pair[0] == pair[1]:
                seeding_type.append("metastatic_confined")
            elif pair[0] != primary and pair[1] == primary:
                seeding_type.append("primary_re_seeding")
            # machina enforces primary tissue to always be source
            # thus any migration from one metastatic tissue to another is part of a cascade
            elif pair[0] != primary and pair[1] != primary and pair[0] != pair[1]:
                #seeding_type.append("cascade_seeding")
                seeding_type.append("metastatic_mono_seeding")
        # check for reseeding or bidirectional seeding
        #aba = find_ABA(asv)
        #if len(aba) > 0:
        #    for pattern in aba:
        #        if pattern[0] == primary:
        #            seeding_type.append("primary_reseeding")
        #        # ignore met->pri->met pattern
        #        elif primary not in pattern:
        #            seeding_type.append("metastatic_reseeding")
        #else:
        #    pass
    return(seeding_type)



# read in big results file for all CP
# output one line per CP as well as one total line
# each line includes summary info as well as each cell in the transition table

parent_dict = {}
leaf_list = []
col_dict = {}
lab_dict = {}
migrations = []
cur_cp = ""
ptissue = sys.argv[2]

tissues = set()
cp = set()

with open(sys.argv[1],'r') as infile:
    lines = infile.readlines()
    for l in lines:
        l = l.strip()
        l = l.split(" ")
        cp.add(l[0])
        if l[1] == "color":
            tissues.add(l[2])

td = {}
for c in cp:
    td[c] = {}
    for t1 in tissues:
        td[c][t1] = {}
        for t2 in tissues:
            td[c][t1][t2] = 0

# print header
#'cascade_seeding': 35, 'metastatic_reseeding': 0, 'metastatic_seeding': 35, 'metastatic_stationary': 73,
#'parallel_seeding': 0, 'primary_reseeding': 0, 'primary_seeding': 72, 'primary_stationary
tabular_tree = []
tissue_dict = {}
print("CP,metastatic_confined,metastatic_mono_seeding,metastatic_parallel_seeding,metastatic_re_seeding,primary_confined,primary_mono_seeding,primary_parallel_seeding,primary_re_seeding")
with open(sys.argv[1],'r') as infile:
    lines = infile.readlines()
    lines.append("END")
    for l in lines:
        l = l.strip()
        l = l.split(" ")
        cp = l[0]
        if cur_cp == "":
            cur_cp = cp
        if cp != cur_cp or l == ["END"]:
            d = td[cur_cp]
            cp_list = []
            topologies = []
            for asv in leaf_list:
                if "XXX" not in asv:
                    parent_list = []
                    cur_parent = list(parent_dict.keys())[get_index(list(parent_dict.values()),asv)] 
                    parent_list.append(cur_parent)
                    while cur_parent != "0" and cur_parent != "0_0" and cur_parent != "grey":
                        cur_parent = list(parent_dict.keys())[get_index(list(parent_dict.values()),cur_parent)] 
                        parent_list.append(cur_parent)
                    parent_tissues = []
                    child_tissue = col_dict[lab_dict[asv]]
                    for parent in reversed(parent_list):
                        parent_tissue = col_dict[lab_dict[parent]]
                        parent_tissues.append(parent_tissue)
                    parent_tissues.append(child_tissue)
                    topologies.append(parent_tissues)
            
            #topo = seeding_topology(topologies,ptissue)
            topo_tree = seeding_topology_tree(tabular_tree,tissue_dict,ptissue)
            # reset tabular tree
            tabular_tree = []
            tissue_dict = {}

            top_types = ["primary_mono_seeding",
                    "primary_re_seeding",
                    "primary_parallel_seeding",
                    "metastatic_parallel_seeding",
                    "metastatic_mono_seeding",
                    "metastatic_re_seeding",
                    "primary_confined",
                    "metastatic_confined"]

            topo_counts = Counter(topo_tree)
            for top_type in top_types:
                if top_type not in topo_counts:
                    topo_counts[top_type] = 0
            topo_counts = dict(sorted(topo_counts.items()))
            ## OUTPUT LINE
            # output values are sorted alphabetically by topology name
            outline = [cur_cp] + list(topo_counts.values())
            print(*outline,sep=",")
            col_dict = {}
            lab_dict = {}
            parent_dict = {}
            leaf_list = []
            model = ""
            migrations = []
            cur_cp = cp
            if l == ["END"]:
                break

        group = l[1]
        distance = 1
        if group == "color":
            col_dict[l[3]] = l[2]
        #elif group == "model":
        #    model = l[2]
        #elif group == "migration":
        #    migrations.append([l[2],l[3]])
        elif group == "label":
            lab_dict[l[3]] = l[2]
        elif group == "tree":
            child = l[3]
            parent = l[2]
            parent_col = lab_dict[parent]
            child_col = lab_dict[child]
            parent_tis = col_dict[parent_col]
            child_tis = col_dict[child_col]

            if parent not in tissue_dict:
                tissue_dict[parent] = parent_tis
            if child not in tissue_dict:
                tissue_dict[child] = child_tis
            if "XXX" not in child:
                tabular_tree.append((parent,child,distance))
            if "ASV" in l[3]:
                leaf = True
                leaf_list.append(l[3])
            else:
                leaf = False
            if parent_tis == child_tis and leaf:
                pass
            else:
                pass
                #td[cp][parent_tis][child_tis] += 1
                #td['all'][parent_tis][child_tis] += 1
      
            if l[2] in parent_dict:
                parent_dict[l[2]].append(l[3])
            else:
                parent_dict[l[2]] = [l[3]]
