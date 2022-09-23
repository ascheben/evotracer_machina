import sys
from collections import Counter
from scipy.stats import entropy

def weird_division(n, d):
    return n / d if d else 0

def get_index(mylist,myname):
    for i,n in enumerate(mylist):
        if myname in n:
            return i

#def seeding_topology(seeding_path):
#    primary = "PRL"
#    # collapse neighbors that are same tissue
#    spath = []
#    cur_tissue = ""
#    for e in seeding_path:
#        if e != cur_tissue:
#            spath.append(e)
#            cur_tissue = e
#    if len(spath) == 2 and spath[0] == primary:
#        seeding_type = "primary_seeding"
#    elif len(spath) == 1 and spath[0] == primary:
#        seeding_type = "no seeding"
#    elif len(spath) == 2 and spath[0] != primary and spath[1] == primary:
#        seeding_type = "re-seeding"
#    elif len(spath) > 2 and spath[0] == primary and primary not in spath[1:]:
#        seeding_type = "seeding cascade"
#    else:
#        seeding_type = "unknown"
#
#    return(seeding_type)

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
        #print("neighbors",neighbors)
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
                #print(f"Adding {spath[-1]} to pri_cas_tissues")
                #print(spath,"seeding cascade")
            # check for reseeding or bidirectional seeding
            aba = find_ABA(asv)
            if len(aba) > 0:
                #print(spath,"found ABA")
                for pattern in aba:
                    if pattern[0] == primary:
                        seeding_type.append("reseeding")
                    else:
                        seeding_type.append("bidirectional_seeding")
            else:
                pass
                #print(spath,"No ABA found")
        #else:
        #    seeding_type = "unknown"

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
model = ""
migrations = []
cur_cp = ""

tissues = ["PRL","LGR","HMR"]

cp = ["CP01","CP02","CP03","CP04","CP05","CP06","CP07","CP08","CP09","CP10","CP11","CP12","CP13","CP14","CP15","CP16","CP17","CP18","CP19","CP20","CP21","CP22","CP23","CP24","CP25","CP26","CP27","CP28","CP29","CP30","CP31","CP32","CP33","CP34","CP35","CP36","CP37","CP38","CP39","CP40","CP41","CP42","CP43","CP44","CP45","CP46","CP47","CP48","CP49","CP50","CP51","CP52","CP53","CP54","CP55","CP56","CP57","CP58","CP59","CP60","CP61","all"]
#cp = ["CP01","CP02","CP03","CP04","CP05","CP06","CP07","CP08","CP09","CP10","CP11","CP12","CP13","CP14","CP15","CP16","CP17","CP18","CP19","CP20","CP21","CP22","CP23","CP24","CP25","CP26","CP27","CP28","CP30","CP31","CP32","CP34","CP35","CP36","CP38","CP39","CP42","CP43","CP49","CP50","CP51","CP52","CP56","CP59","all"]
td = {}
for c in cp:
    td[c] = {}
    for t1 in tissues:
        td[c][t1] = {}
        for t2 in tissues:
            td[c][t1][t2] = 0

with open(sys.argv[1],'r') as infile:
    for l in infile:
        l = l.strip()
        l = l.split(" ")
        cp = l[0]
        print("CP position", cp, cur_cp)
        if cur_cp == "":
            cur_cp = cp
        if cp != cur_cp:
            # output CP summary
            #print("Finished CP",cur_cp)
            # calculate entropy
            d = td[cur_cp]
            source_matrix = []
            source_entropy = []
            all_tis_source_tot = 0
            all_source_probs = []
            for t1 in tissues:
                source_tot = 0
                target_tot = 0
                source_probs = []
                target_probs = []
                for t2 in tissues:
                    source_tot += d[t1][t2]
                    target_tot += d[t2][t1]
                    all_tis_source_tot += d[t1][t2]
                for t2 in tissues:
                     source_probs.append(weird_division(d[t1][t2],source_tot))
                     target_probs.append(weird_division(d[t2][t1],target_tot))
                hsource = entropy(source_probs,base=3)
                htarget = entropy(target_probs,base=3)
                source_matrix.append(source_probs)
                source_entropy.append(hsource)
            for t1 in tissues:
                for t2 in tissues:
                    all_source_probs.append(weird_division(d[t1][t2],all_tis_source_tot))

            #print(cur_cp,model,len(migrations),source_entropy[0],source_entropy[1],source_entropy[2],source_matrix,all_source_probs)

            
            print("leaf list:", leaf_list)
            #print("parent dict:",parent_dict)
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
                    #topo = seeding_topology(parent_tissues)
                    topologies.append(parent_tissues)
                    print(cur_cp,asv,list(reversed(parent_list)),parent_tissues)
                #cp_list.append(topo)
            topo = seeding_topology(topologies)
            #print("topologies,topo",topologies,topo)
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
                    
            print(cur_cp, col_dict,topo_counts)
            col_dict = {}
            lab_dict = {}
            parent_dict = {}
            leaf_list = []
            model = ""
            migrations = []
            cur_cp = cp
        group = l[1]
        #print(l)
        if group == "color":
            col_dict[l[3]] = l[2]
        elif group == "model":
            model = l[2]
        elif group == "migration":
            migrations.append([l[2],l[3]])
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
            #print("l,parent_tis,child_tis",l,parent_tis,child_tis,leaf)
            if parent_tis == child_tis and leaf:
                #print("no transition")
                pass
            #if cp not in trans_dict:
            #    trans_dict[cp] = {}
            #if parent_tis not in trans_dict[cp]:
            #    trans_dict[cp][parent_tis] = {}
            #if child_tis not in trans_dict[cp][parent_tis]:
            #    trans_dict[cp][parent_tis][child_tis] = 0
            else:
                #print("Counting transition from ",parent_tis," to ",child_tis)
                td[cp][parent_tis][child_tis] += 1
                td['all'][parent_tis][child_tis] += 1
      
            if l[2] in parent_dict:
                parent_dict[l[2]].append(l[3])
            else:
                parent_dict[l[2]] = [l[3]]
                
# print entropy for all
d = td['all']
source_matrix = []
source_entropy = []
all_tis_source_tot = 0
all_source_probs = []

for t1 in tissues:
    source_tot = 0
    target_tot = 0
    source_probs = []
    target_probs = []
    for t2 in tissues:
        source_tot += d[t1][t2]
        target_tot += d[t2][t1]
        all_tis_source_tot += d[t1][t2]
    for t2 in tissues:
         source_probs.append(weird_division(d[t1][t2],source_tot))
         target_probs.append(weird_division(d[t2][t1],target_tot))
  
    hsource = entropy(source_probs,base=3)
    htarget = entropy(target_probs,base=3)
    source_matrix.append(source_probs)
    source_entropy.append(hsource)
for t1 in tissues:
    for t2 in tissues:
        all_source_probs.append(weird_division(d[t1][t2],all_tis_source_tot))

#print("cols:",col_dict)
#print("labels:", lab_dict)
#print("td:",td)
#print("all","NA","NA",source_entropy,source_matrix,all_source_probs)

