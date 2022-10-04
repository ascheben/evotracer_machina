import sys
from ete3 import Tree
import math
from math import comb
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
    if fN < N and fk < k and fN > 0 and fN > fk and N > k:
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

def traverse_tree(tabular_tree,freq_dict,cp):
    #print(tabular_tree)
    tree = Tree.from_parent_child_table(tabular_tree)
    #print(tree.write())
    # ignore root node
    for node in tree.iter_descendants("levelorder"):
        if not node.is_leaf():
            sisters = []
            sister_descendants = []
            sister_descendant_freqs = []
            for n in node.get_sisters():
                sisters.append(n.name)
                descendants = n.get_leaf_names()
                sister_descendants = sister_descendants + descendants
                for d in descendants:
                    if d in freq_dict[cp]:
                        sister_descendant_freqs.append(freq_dict[cp][d])
            node_descendants = []
            node_descendant_freqs = []
            for nd in node.get_leaf_names():
                node_descendants.append(nd)
                if nd in freq_dict[cp]:
                    node_descendant_freqs.append(freq_dict[cp][nd])
            #print(node.name,sisters,leaf_descendants,sister_descendant_freqs)
    
            k = 1 + len(sisters)
            fk = 1
            N = len(node_descendants) + len(sister_descendants)
            fn = len(node_descendants)
            N_asv = sum(node_descendant_freqs) + sum(sister_descendant_freqs)
            fn_asv = sum(node_descendant_freqs)
            pval = log_pvalue(k, fk,N,fn)
            pval_asv = log_pvalue(k, fk,N_asv,fn_asv)
            print(cp,node.name,k,fk,N,fn,N_asv,fn_asv,pval,pval_asv)

# get asv freq information from file
# this may be needed if no information is availble in the tree file
if len(sys.argv) == 3:
    asvdict = {}
    with open(sys.argv[2],'r') as asv:
        header = next(asv)
        for l in asv:
            l = l.strip()
            l = l.split(",")
            name = l[0]
            tissue = l[1]
            count = int(l[3])
            cp_name = l[29]
            name_tissue = name + "_" + tissue
            if cp_name not in asvdict:
                asvdict[cp_name] = {}
            asvdict[cp_name][name_tissue] = count
            # for single tissue ASVs
            if name not in asvdict[cp_name]:
                asvdict[cp_name][name] = count

tabular_tree = []
freq_dict = {}
header = ["CP","clade","total_lineages","clade_lineages","total_leaves","clade_leaves","total_asv_freq","clade_asv_freq","relate_logp","relate_logp_asv_freq"]
print(*header,sep=" ")
cur_cp = ""
with open(sys.argv[1],'r') as infile:
    lines = infile.readlines()
    lines.append("END")
    for l in lines:
        l = l.strip()
        l = l.split(" ")
        if l[0] == "END":
            traverse_tree(tabular_tree,freq_dict,cp)
        elif l[1] == "tree":
            cp = l[0]
            if cur_cp == "":
                cur_cp = cp
            if cp != cur_cp:
                traverse_tree(tabular_tree,freq_dict,cur_cp)
                tabular_tree = []
                cur_cp = cp
            elif cur_cp == cp:
                parent = str(l[2])
                child = str(l[3])
                distance = 1
                try:
                    freq = l[4]
                except:
                    if child in asvdict[cp]:
                        freq = asvdict[cp][child]
                    else:
                        freq = "NA"
                if cp not in freq_dict:
                    freq_dict[cp] = {}
                if freq != "NA":
                    freq_dict[cp][child] = int(freq)
                if "XXX" not in child:
                    tabular_tree.append((parent,child,distance))
            #elif l[1] != "tree" and len(tabular_tree) > 0:
            #    traverse_tree(tabular_tree,freq_dict,cp)
            #    tabular_tree = []
# expected input file
#CP02 tree 0 0^LVM NA 0.36364
#CP02 tree 0 0^RBL NA 0.27273
#CP02 tree 0 ASV508 96 0.36364

