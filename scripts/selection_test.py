import sys
import math
from math import comb
from collections import defaultdict
from collections import Counter
import numpy as np

def weird_division(n, d):
    return n / d if d else 0

def get_index(mylist,myname):
    for i,n in enumerate(mylist):
        if myname in n:
            return i

def log_pvalue(k, fk,N,fN):
    '''
    k = all lineages
    fk = all target lineages
    N = all leaves
    fn = all target leaves
    '''
    if fN < N and fk < k and fN > 0:
        px  = math.log(math.factorial(N-fN-1)) - math.log(math.factorial(k-fk-1)) - math.log(math.factorial(N-k+fk-fN))
        px += math.log(math.factorial(fN-1)) - math.log(math.factorial(fk-1)) - math.log(math.factorial(fN-fk))
        px -= math.log(math.factorial(N-1)) - math.log(math.factorial(k-1)) - math.log(math.factorial(N-k))
        logp = px
        x = fN - fk
        y = N - k
        c = N - 1
        while x < N-k:
            var  = fk + x
            px  += np.log( (y-x)/(x+1.0) * var / (c - var) )
            logp = np.log( 1.0 + np.exp( px - logp ) ) + logp;
            x += 1
        if logp > 0:
            logp = 0
        return(logp)
    else:
        return("NA")

def yang_pvalue(k, fk,N,fN):
    #nc = fN
    sum_p_Nk_nc = 0
    for nc in range(fn,N-k):
        p_Nk_nc = math.comb(N-nc-1,k-2) / math.comb(N-1,k-1)
        sum_p_Nk_nc += p_Nk_nc
    if sum_p_Nk_nc == 0:
        return -np.inf
    else:
        return math.log(sum_p_Nk_nc)



tissues = set()
cp_set = set(["all"])
asv_dict = defaultdict(lambda: defaultdict(int))
with open(sys.argv[2],'r') as infile:
    header = next(infile)
    for l in infile:
        l = l.strip()
        l = l.split(",")
        asv = l[0]
        tissue = l[1]
        tissues.add(tissue)
        mutation = l[25]
        cp_name = l[29]
        cp_set.add(cp_name)
        asv_count = int(l[3])
        asv_dict[cp_name][asv+"_"+tissue] = asv_count


parent_dict = {}
leaf_list = []
col_dict = {}
lab_dict = {}
model = ""
migrations = []
cur_cp = ""


td = {}
for c in cp_set:
    td[c] = {}
    for t1 in tissues:
        td[c][t1] = {}
        for t2 in tissues:
            td[c][t1][t2] = 0
# print header
print("CP,clade,clade_leaves,sister_leaves,other_leaves,clade_samples,sister_samples,nonclade_samples,clade_lineages,sister_lineages,nonclade_lineages,sister_relate_logp,nonclade_relate_logp,sister_sample_relate_logp,nonclade_sample_relate_logp")
with open(sys.argv[1],'r') as infile:
    lines = infile.readlines()
    lines.append("END")
    for l in lines:
        l = l.strip()
        l = l.split(" ")
        cp = l[0]
        if cur_cp == "":
            cur_cp = cp
        if cp != cur_cp or cp == ["END"]:
            # output CP summary
            # calculate entropy
            d = td[cur_cp]

            cp_list = []
            topologies = []
            for asv in leaf_list:
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
                #print(cur_cp,asv,list(reversed(parent_list)),parent_tissues)
                #cp_list.append(topo)

            count_dict = {}
            for k,v in parent_dict.items():
                for asv in v:
                    if "XXX" in asv:
                        count_dict["ASVXXX"] = 0
                    else:
                        try:
                            if asv_dict[cur_cp][asv] != 0:
                                #print(asv,asv_dict[cur_cp][asv])
                                count_dict[asv] = asv_dict[cur_cp][asv]
                            else:
                                for suffix in tissues:
                                    try:
                                        if asv_dict[cur_cp][asv+"_"+suffix] != 0:
                                            #print(asv,asv_dict[cur_cp][asv+"_"+suffix])
                                            count_dict[asv+"_"+suffix] = asv_dict[cur_cp][asv+"_"+suffix]
                                    except:
                                        pass
                        except:
                            pass
            #print(cur_cp,parent_dict,count_dict)
            # Test for selection
            child_dict = defaultdict(set)
            for k,v in parent_dict.items():
                lowest_node = True
                for node in v:
                    if not "ASV" in node:
                        lowest_node = False
                if lowest_node:
                    traversed = []
                    cur_node = k
                    while cur_node != "0":
                        cur_parent = list(parent_dict.keys())[get_index(list(parent_dict.values()),cur_node)]
                        if cur_parent != "0":
                            child_dict[cur_node].add(cur_parent)
                        for t in traversed:
                            if t != "0":
                                if cur_parent != "0":
                                    child_dict[t].add(cur_parent)
                        traversed.append(cur_node)
                        cur_node = cur_parent

            for key in parent_dict:
                if key != "0":
                    sisters = set()
                    # find parent of key
                    key_parent = list(parent_dict.keys())[get_index(list(parent_dict.values()),key)]
                    # find all children of parent
                    for node in parent_dict[key_parent]:
                        if node != key and "ASV" not in node and node != "0" and node != key_parent:
                            sisters.add(node)
                            child_nodes = child_dict[node]
                            for c in child_nodes:
                                if c != key and "ASV" not in c and c != "0" and c != key_parent:
                                    sisters.add(c)
                    #print(key, " parent is ", key_parent, " with sisters: ",sisters)
                    total_leaves = 0
                    total_leaf_samples = 0 
                    for value in parent_dict[key]:
                        if "ASV" in value and "ASVXXX" not in value:
                            total_leaf_samples += 1
                            if value in count_dict:
                                total_leaves += count_dict[value]
                            else:
                                for suffix in tissues:
                                    try:
                                        nleaves = count_dict[value+"_"+suffix]
                                    except:
                                        pass
                                total_leaves += nleaves
                    total_leaves_sister = 0
                    total_leaf_samples_sister = 0
                    for s in sisters:
                        for value in parent_dict[s]:
                            if "ASV" in value and "ASVXXX" not in value:
                                total_leaf_samples_sister +=1
                                if value in count_dict:
                                    total_leaves_sister += count_dict[value]
                                else:
                                    for suffix in tissues:
                                        try:
                                            nleaves = count_dict[value+"_"+suffix]
                                        except:
                                            pass
                                    total_leaves_sister += nleaves
                    

                    total_leaves_other = 0
                    total_leaf_samples_other = 0
                    total_other = 0
                    for s in parent_dict:
                        if s != "0" and s != key :
                            total_other += 1
                            for value in parent_dict[s]:
                                if "ASV" in value and "ASVXXX" not in value:
                                    total_leaf_samples_other += 1
                                    if value in count_dict:
                                        total_leaves_other += count_dict[value]
                                    else:
                                        for suffix in tissues:
                                            try:
                                                nleaves = count_dict[value+"_"+suffix]
                                            except:
                                                pass
                                        total_leaves_other += nleaves
     

                    #print("len(sisters),total_leaves,total_leaves_sister,total_other,total_leaves_other:",len(sisters),total_leaves,total_leaves_sister,total_other,total_leaves_other)         

                    # Use only lineages with same parent as target
                    k = len(sisters) + 1
                    fk = 1
                    N = total_leaves + total_leaves_sister
                    fn = total_leaves

                    sis_relate_logp = log_pvalue(k,fk,N,fn)
                    sis_yang_logp = yang_pvalue(k,fk,N,fn)

                    # Use whole tree
                    k = total_other + 1
                    fk = 1
                    N = total_leaves + total_leaves_other
                    fn = total_leaves

                    all_relate_logp = log_pvalue(k,fk,N,fn)
                    all_yang_logp = yang_pvalue(k,fk,N,fn)

                    # Use lineage counts not asv counts for leaves
                    k = len(sisters) + 1
                    fk = 1
                    N = total_leaf_samples + total_leaf_samples_sister
                    fn = total_leaf_samples

                    sis_samples_relate_logp = log_pvalue(k,fk,N,fn)

                    k = total_other + 1
                    fk = 1
                    N = total_leaf_samples + total_leaf_samples_other 
                    fn = total_leaf_samples

                    all_samples_relate_logp = log_pvalue(k,fk,N,fn)
 
                    print(*[cur_cp,key,total_leaves,total_leaves_sister,total_leaves_other,total_leaf_samples,total_leaf_samples_sister,total_leaf_samples_other,1,len(sisters),total_other,sis_relate_logp,all_relate_logp,sis_samples_relate_logp,all_samples_relate_logp],sep=",")
                    

            # reset data
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
                
