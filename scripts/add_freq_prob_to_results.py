import sys
from collections import defaultdict

pdict = {}
with open(sys.argv[1],'r') as mat:
    header = next(mat)
    for l in mat:
        l = l.strip()
        l = l.split(" ")
        pdict[l[0]] = {}
        prob_mat = [l[6:9],l[9:12],l[12:]]
        tis_list = ["PRL","LGR","HMR"]
        for i,tis in enumerate(tis_list):
            pdict[l[0]][tis] = {}
            for k,tis2 in enumerate(tis_list):
                trans_prob = prob_mat[i][k]
                pdict[l[0]][tis][tis2] = trans_prob
        #print(l[0],pdict[l[0]])
print(pdict)
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

cur_cp = ""
coldict = {}
labdict = {}
with open(sys.argv[3],'r') as results:
    header = next(results)
    for l in results:
        l = l.strip()
        l = l.split(" ")
        if l[0] != cur_cp:
            cur_cp = l[0]
        if l[1] == "color":
            if cur_cp not in coldict:
                coldict[cur_cp] = {}
            coldict[cur_cp][l[3]] = l[2]
        if l[1] == "label":
            if cur_cp not in labdict:
                labdict[cur_cp] = {}
            labdict[cur_cp][l[3]] = l[2]

        if l[1] == "tree":
            parent = l[2]
            child = l[3]

            source = coldict[cur_cp][labdict[cur_cp][parent]]
            target = coldict[cur_cp][labdict[cur_cp][child]]

            if "ASV" in child and not "XXX" in child:
                child_freq = asvdict[cur_cp][child]
            else:
                child_freq = "NA"
            trans_prob = pdict[cur_cp][source][target]
            #print(cur_cp,parent,child,source,target,child_freq,trans_prob)

with open(sys.argv[3],'r') as results:
    header = next(results)
    for l in results:
        l = l.strip()
        l = l.split(" ")
        if l[0] != cur_cp:
            cur_cp = l[0]
        if l[1] == "tree":
            parent = l[2]
            child = l[3]
            source = coldict[cur_cp][labdict[cur_cp][parent]]
            target = coldict[cur_cp][labdict[cur_cp][child]]

            if "ASV" in child and not "XXX" in child:
                child_freq = asvdict[cur_cp][child]
            else:
                child_freq = "NA"
            trans_prob = pdict[cur_cp][source][target]
            outline = l + [child_freq,trans_prob]
            print(*outline,sep=" ")
        else:
            outline = l + ["NA","NA"]
            print(*outline,sep=" ")

