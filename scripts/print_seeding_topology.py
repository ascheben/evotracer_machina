import sys
from collections import Counter

def get_index(mylist,myname):
    for i,n in enumerate(mylist):
        if myname in n:
            return i

def find_ABA(seeding_path):
    '''
    find patterns of tissue transitions indicating
    bidirectional seeding or reseeding
    '''
    ABA = []
    for i,s in enumerate(seeding_path):
        try:
            second = seeding_path[i+1]
            third = seeding_path[i+2]
            if i == third and i != second:
                pattern = [i,second,third]
                ABA.append(pattern)
        except:
            pass
    return ABA


def seeding_topology(seeding_path_list):
    '''
    Takes in nested list of seeding paths for each ASV
    '''
    primary = "PRL"
    # asv tissue for primary and cascade
    pri_cas_tissues = set()
    pri2met_tissues = set()
    seeding_type = []
    for asv in seeding_path_list: 
        #spath = []
        # collapse neighbors that are same tissue
        #cur_tissue = ""
        #for e in asv:
        #    if e != cur_tissue:
        #        spath.append(e)
        #        cur_tissue = e
        # seeding is always from primary in this case
        #if len(spath) == 2 and spath[0] == primary:
        #    seeding_type.append("primary_seeding")
        #    pri_cas_tissues.add(spath[1])
        #    pri2met_tissues.add(spath[1])
        #elif len(spath) == 1 and spath[0] == primary:
        #    seeding_type.append("primary_expansion")
        #elif len(spath) > 2:
        #    if spath[0] == primary and primary not in spath[1:]:
        #seeding_type.append("seeding_cascade")
        #pri_cas_tissues.add(spath[-1])

        window = 2
        overlap = 1
        neighbors = [asv[i:i+window] for i in range(0, len(asv), window-overlap)]
        for pair in neighbors:
            if len(pair) == 2:
                if pair[0] == primary and pair[1] != primary:
                    pri2met_tissues.add(pair[1])
                    seeding_type.append("primary_seeding")
                elif pair[0] == primary and pair[1] == primary:
                    seeding_type.append("primary_expansion")
                elif pair[0] != primary and pair[1] != primary and pair[0] == pair[1]:
                    seeding_type.append("metastatic_expansion")
                # machina enforces primary tissue to always be source
                # thus any migration from one metastatic tissue to another is part of a cascade
                elif pair[0] != primary and pair[1] != primary and pair[0] != pair[1]:
                    seeding_type.append("seeding_cascade")
            # check for reseeding or bidirectional seeding
            aba = find_ABA(asv)
            if len(aba) > 0:
                for pattern in aba:
                    if pattern[0] == primary:
                        seeding_type.append("reseeding")
                    else:
                        seeding_type.append("bidirectional_seeding")
            else:
                pass
    # check for parallel seeding
    #if len(pri_cas_tissues) > 1:
    if len(pri2met_tissues) > 1:
        seeding_type.append("parallel_seeding")

    return(seeding_type)




#entropy([5/11, 6/11, 0/11])

# read in big results file for all CP
# output one line per CP as well as one total line
# each line includes summary info as well as each cell in the transition table

parent_dict = {}
leaf_list = []
col_dict = {}
lab_dict = {}
migrations = []
cur_cp = ""


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
                    while cur_parent != "0":
                        cur_parent = list(parent_dict.keys())[get_index(list(parent_dict.values()),cur_parent)] 
                        parent_list.append(cur_parent)
                    parent_tissues = []
                    child_tissue = col_dict[lab_dict[asv]]
                    for parent in reversed(parent_list):
                        parent_tissue = col_dict[lab_dict[parent]]
                        parent_tissues.append(parent_tissue)
                    parent_tissues.append(child_tissue)
                    topologies.append(parent_tissues)
                    #print(cur_cp,asv,list(reversed(parent_list)),parent_tissues)
            topo = seeding_topology(topologies)
            top_types = ["seeding_cascade",
                    "primary_seeding",
                    "reseeding",
                    "bidirectional_seeding",
                    "parallel_seeding",
                    "primary_expansion",
                    "metastatic_expansion"]
            topo_counts = Counter(topo)
            for top_type in top_types:
                if top_type not in topo_counts:
                    topo_counts[top_type] = 0
            topo_counts = dict(sorted(topo_counts.items()))
                    
            #print(cur_cp, col_dict,topo_counts)
            ## OUTPUT LINE
            #outline = [cur_cp] + list(topo_counts.values()) + list(topo_counts.keys())
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
        if group == "color":
            col_dict[l[3]] = l[2]
        #elif group == "model":
        #    model = l[2]
        #elif group == "migration":
        #    migrations.append([l[2],l[3]])
        elif group == "label":
            lab_dict[l[3]] = l[2]
        elif group == "tree":
            parent_col = lab_dict[l[2]]
            child_col = lab_dict[l[3]]
            parent_tis = col_dict[parent_col]
            child_tis = col_dict[child_col]
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
                

