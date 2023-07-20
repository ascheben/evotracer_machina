import sys
from collections import defaultdict

pdict = {}
with open(sys.argv[1],'r') as mat:
    header = next(mat)
    header = header.strip()
    header = header.split(",")
    # find number of tissues
    #CP,model,migrations,edges,TreeMetRate,H_HMR,H_LGR,H_PRL,HMR_HMR,HMR_LGR,HMR_PRL,LGR_HMR,LGR_LGR,LGR_PRL,PRL_HMR,PRL_LGR,PRL_PRL
    # CP,model,migrations,edges,TreeMetRate,H_LVM,H_LVR,H_PRL,H_RBL,LVM:LVM,LVM:LVR,LVM:PRL,LVM:RBL,LVR:LVM,LVR:LVR,LVR:PRL,LVR:RBL,PRL:LVM,PRL:LVR,PRL:PRL,PRL:RBL,RBL:LVM,RBL:LVR,RBL:PRL,RBL:RBL
    tissue_set = set()
    for elem in header[5:]:
        if elem.startswith("H_"):
            tis_name = elem.replace("H_","")
            tissue_set.add(tis_name)
    tissue_cnt = len(tissue_set)
    prob_index = 5 + tissue_cnt
    for l in mat:
        l = l.strip()
        l = l.split(",")
        pdict[l[0]] = {}
        #prob_mat = [l[8:11],l[11:14],l[14:]]
        # script breaks between MMUS1495 and 1469 
        prob_mat = l[prob_index:]
        #prob_mat = l[9:]
        #tis_list = ["PRL","LGR","HMR"]
        prob_names = header[prob_index:]
        #prob_names = header[9:]
        for i,tis in enumerate(prob_names):
            tis1 = tis.split(":")[0]
            tis2 = tis.split(":")[1]
            if tis1 not in pdict[l[0]]:
                pdict[l[0]][tis1] = {}
            if tis2 not in pdict[l[0]][tis1]:
                pdict[l[0]][tis1][tis2] = {}
            pdict[l[0]][tis1][tis2] = int(prob_mat[i])

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
    #header = next(results)
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
    #header = next(results)
    header = "CP group value1 value2 asv_freq_count transition_probability transition_probability_no_selfmigration"
    print(header)
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
            trans_count = pdict[cur_cp][source][target]
            total_edges = 0
            total_nonself_edges = 0
            for k1,v1 in pdict[cur_cp].items():
                for k2,v2 in v1.items():
                    total_edges += v2
                    if k1 != k2:
                        total_nonself_edges += v2

            trans_prob_all = round(float(trans_count / total_edges),3)
            if source != target:
                trans_prob_nonself = round(float(trans_count / total_nonself_edges),3)
            else:
                trans_prob_nonself = 0
            outline = l + [child_freq,trans_prob_all,trans_prob_nonself]
            print(*outline,sep=" ")
        else:
            outline = l + ["NA","NA","NA"]
            print(*outline,sep=" ")

