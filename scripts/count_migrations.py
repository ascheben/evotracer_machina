import sys
from scipy.stats import entropy

def weird_division(n, d):
    return n / d if d else 0

#entropy([5/11, 6/11, 0/11])

# read in big results file for all CP
# output one line per CP as well as one total line
# each line includes summary info as well as each cell in the transition table


col_dict = {}
lab_dict = {}
model = ""
migrations = []
cur_cp = ""

tissues = ["PRL","LGR","HMR"]
cp = ["CP01","CP02","CP03","CP04","CP05","CP06","CP07","CP08","CP09","CP10","CP11","CP12","CP13","CP14","CP15","CP16","CP17","CP18","CP19","CP20","CP21","CP22","CP23","CP24","CP25","CP26","CP27","CP28","CP30","CP31","CP32","CP34","CP35","CP36","CP38","CP39","CP42","CP43","CP49","CP50","CP51","CP52","CP56","CP59","all"]
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

            print(cur_cp,model,len(migrations),source_entropy[0],source_entropy[1],source_entropy[2],source_matrix,all_source_probs)
            col_dict = {}
            lab_dict = {}
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



print("all","NA","NA",source_entropy,source_matrix,all_source_probs)

