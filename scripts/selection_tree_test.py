#!/usr/bin/env python3

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



def log_pvalue(k, fk,N,fN):
    '''
    k = all lineages
    fk = all target lineages
    N = all leaves
    fn = all target leaves
    '''
    if max(k, fk,N,fN) > 100000:
        return "NA"
    elif fN < N and fk < k and fN > 0 and fN > fk and N > k:
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
    sig_threshold = math.log(0.1)
    odds_threshold = 2
    tree = Tree.from_parent_child_table(tabular_tree)
    #print(tree.write())
   
    tree_depth = get_tree_depth(tree)

    expanded_asv_lineages = set()
    expanded_asv_counts = set()
    expanded_asv_lineages_relax = set()
    expanded_asv_counts_relax = set()
       
    all_asv = set()
    nodes_tested = set()
    nodes_tested_relax = set()
    nodes_total = set()
    nodes_expanded = 0
    nodes_expanded_asv = 0
    nodes_expanded_relax = 0
    nodes_expanded_asv_relax = 0

    # fully unresolved trees should still produce one "test" line
    if tree_depth == 1:
        print(*["test",cp,"NA","NA","NA","NA","NA","NA","NA","NA","NA"],sep="\t")
    # ignore root node
    for node in tree.iter_descendants("levelorder"):
        if not node.is_leaf():
            nodes_total.add(node.name)
            sisters = []
            sister_descendants = []
            sister_descendant_freqs = []
            # allow single tip sisters
            sisters_relax = []
            sister_descendants_relax = []
            sister_descendant_freqs_relax = []

            for n in node.get_sisters():
                sisters_relax.append(n.name)
                descendants_relax = n.get_leaf_names()
                sister_descendants_relax = sister_descendants_relax + descendants_relax
                for d in descendants_relax:
                    if d in freq_dict[cp]:
                        sister_descendant_freqs_relax.append(freq_dict[cp][d])
                # ignore sister lineages that have just one leaf
                if not n.is_leaf():
                    sisters.append(n.name)
                    descendants = n.get_leaf_names()
                    sister_descendants = sister_descendants + descendants
                    for d in descendants:
                        if d in freq_dict[cp]:
                            sister_descendant_freqs.append(freq_dict[cp][d])
            if len(sisters)>0:
                nodes_tested.add(node.name)
            if len(sisters_relax)>0:
                nodes_tested_relax.add(node.name)


            node_descendants = []
            node_descendant_freqs = []
            for nd in node.get_leaf_names():
                node_descendants.append(nd)
                if nd in freq_dict[cp]:
                    node_descendant_freqs.append(freq_dict[cp][nd])
            #print(node.name,sisters,leaf_descendants,sister_descendant_freqs)
            node_descendants_uniq = set()
            for nd in node_descendants:
                head, sep, tail = nd.partition('_')
                node_descendants_uniq.add(head)
            sister_descendants_uniq = set()
            for sd in sister_descendants:
                head, sep, tail = sd.partition('_')
                sister_descendants_uniq.add(head)
            sister_descendants_relax_uniq = set()
            for sdr in sister_descendants_relax:
                head, sep, tail = sdr.partition('_')
                sister_descendants_relax_uniq.add(head)
       
            sisters = list(set([s.partition('_')[0] for s in sisters]))
            sisters_relax = list(set([s.partition('_')[0] for s in sisters_relax]))

            k = 1 + len(sisters)
            fk = 1
            N = len(node_descendants_uniq) + len(sister_descendants_uniq)
            fn = len(node_descendants_uniq)
            N_asv = sum(node_descendant_freqs) + sum(sister_descendant_freqs)
            fn_asv = sum(node_descendant_freqs)
            
            k_relax = 1 + len(sisters_relax)
            N_relax = len(node_descendants_uniq) + len(sister_descendants_relax_uniq)
            N_asv_relax = sum(node_descendant_freqs) + sum(sister_descendant_freqs_relax)

            # Odds ratio for relaxed lineage counts
            oddsratio = ((k_relax-fn)*fn)/(fk*(N_relax-fn))
            #print("k, fk,N,fn",k, fk,N,fn)
            pval = log_pvalue(k, fk,N,fn)
            pval_asv = log_pvalue(k, fk,N_asv,fn_asv)
            pval_relax = log_pvalue(k_relax, fk,N_relax,fn)
            pval_asv_relax = log_pvalue(k_relax, fk,N_asv_relax,fn_asv)
            #if pval != "NA":
            #    if pval < sig_threshold:
            #        expanded_asv_lineages.update(node_descendants)
            #        nodes_expanded += 1
            #if pval_asv != "NA":
            #    if pval_asv < sig_threshold:
            #        expanded_asv_counts.update(node_descendants)
            #        nodes_expanded_asv += 1
            #if pval_relax != "NA":
            #    if pval_relax < sig_threshold:
            #        expanded_asv_lineages_relax.update(node_descendants)
            #        nodes_expanded_relax += 1
            #if pval_asv_relax != "NA":
            #    if pval_asv_relax < sig_threshold:
            #        expanded_asv_counts_relax.update(node_descendants)
            #        nodes_expanded_asv_relax += 1

            if oddsratio > odds_threshold:
                expanded_asv_lineages_relax.update(node_descendants)
                nodes_expanded_relax += 1



            # Print per lineage results
            if len(node_descendants_uniq)>0:
                node_descendants_uniq= ";".join(node_descendants_uniq)
            else: 
                node_descendants_uniq = "NA"
            if len(sister_descendants_uniq)>0:
                sister_descendants_uniq= ";".join(sister_descendants_uniq)
            else: 
                sister_descendants_uniq = "NA"
            if len(sister_descendants_relax_uniq)>0:
                sister_descendants_relax_uniq= ";".join(sister_descendants_relax_uniq)
            else: 
                sister_descendants_relax_uniq = "NA"
            print(*["test",cp,node.name,k_relax,fk,N_relax,fn,oddsratio,pval_relax,node_descendants_uniq,sister_descendants_relax_uniq],sep="\t")
            #print(*["test",cp,node.name,k,fk,N,fn,N_asv,fn_asv,k_relax,N_relax,N_asv_relax,pval,pval_asv,pval_relax,pval_asv_relax,node_descendants_uniq,sister_descendants_uniq,sister_descendants_relax_uniq],sep="\t")

        else:
            leaf_name = node.name
            all_asv.add(leaf_name)
    total_freq = 0
    all_asv_uniq = set()
    for a in all_asv:
        total_freq += freq_dict[cp][a]
        head, sep, tail = a.partition('_')
        all_asv_uniq.add(head)
    eal_count = 0
    if len(expanded_asv_lineages) > 0:
        for a in expanded_asv_lineages:
            eal_count += freq_dict[cp][a]
    eac_count = 0
    if len(expanded_asv_counts) > 0:
        for a in expanded_asv_counts:
            eac_count += freq_dict[cp][a]
    ealr_count = 0
    if len(expanded_asv_lineages_relax) > 0:
        for a in expanded_asv_lineages_relax:
            ealr_count += freq_dict[cp][a]
    eacr_count = 0
    if len(expanded_asv_counts_relax) > 0:
        for a in expanded_asv_counts_relax:
            eacr_count += freq_dict[cp][a]

    # uniquify lists with tissue suffixes
    expanded_asv_lineages = list(set([s.partition('_')[0] for s in expanded_asv_lineages]))
    expanded_asv_counts = list(set([s.partition('_')[0] for s in expanded_asv_counts]))
    expanded_asv_lineages_relax = list(set([s.partition('_')[0] for s in expanded_asv_lineages_relax]))
    expanded_asv_counts_relax = list(set([s.partition('_')[0] for s in expanded_asv_counts_relax]))

    print(*["expansion",cp,tree_depth,len(nodes_total),len(nodes_tested_relax),nodes_expanded_relax,len(all_asv_uniq),len(expanded_asv_lineages_relax),total_freq,ealr_count],sep="\t")    
    #print(*["expansion",cp,len(nodes_total),len(nodes_tested),len(nodes_tested_relax),nodes_expanded,nodes_expanded_asv,nodes_expanded_relax,nodes_expanded_asv_relax,len(all_asv_uniq),len(expanded_asv_lineages),len(expanded_asv_counts),len(expanded_asv_lineages_relax),len(expanded_asv_counts_relax),total_freq,eal_count,eac_count,ealr_count,eacr_count],sep="\t")    
# get asv freq information from file
# this may be needed if no information is availble in the tree file
done_list = []
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
            sequence = l[11]

            done_line = cp_name + name + tissue + sequence
            if done_line not in done_list:
                if cp_name not in asvdict:
                    asvdict[cp_name] = {}
                if name_tissue in asvdict[cp_name]:
                    asvdict[cp_name][name_tissue] = count + asvdict[cp_name][name_tissue]
                else:
                    asvdict[cp_name][name_tissue] = count
                # for single tissue ASVs
                if name not in asvdict[cp_name]:
                    asvdict[cp_name][name] = count
                else:
                    asvdict[cp_name][name] = asvdict[cp_name][name] + count
                done_list.append(done_line)
tabular_tree = []
freq_dict = {}
header_test = ["test","CP","clade","total_lineages","clade_lineages","total_leaves","clade_leaves","odds_ratio","relate_logp","clade_leaf_names","sister_leaf_names"]
#header_test = ["test","CP","clade","total_lineages","clade_lineages","total_leaves","clade_leaves","total_asv_freq","clade_asv_freq", "total_lineages_relax","total_leaves_relax","total_asv_freq_relax","relate_logp","relate_logp_asv_freq","relate_logp_relax","relate_logp_asv_freq_relax","clade_leaf_names","sister_leaf_names","sister_leaf_names_relax"]
print(*header_test,sep="\t")
header_expansion=["expansion","CP","tree_depth","total_nodes","tested_nodes","expanded_nodes","total_leaves","leaves_expanded_asv_lineages","total_asv_freq","expanded_asv_freq"]
#header_expansion=["expansion","CP","total_nodes","tested_nodes","tested_nodes_relax","expanded_nodes","expanded_nodes_asv_freq","expanded_nodes_relax","expanded_nodes_asv_freq_relax","total_leaves","leaves_expanded_asv_lineages","leaves_expanded_asv_counts","leaves_expanded_asv_lineages_relax","leaves_expanded_asv_counts_relax","total_asv_freq","eal_count","eac_count","eal_relax_count","eal_relax_count"]
print(*header_expansion,sep="\t")


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

