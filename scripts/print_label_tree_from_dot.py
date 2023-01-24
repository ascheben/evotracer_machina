import sys
import re


#grep "\->" T-PRL-S-binarized.dot| cut -d' ' -f1,3|sed 's/^\t//' > T-PRL-S-binarized_parent_child.txt
#grep label T-PRL-S-binarized.dot| tr ',' ' '| cut -d' ' -f1,4,5| sed 's/^\t//'| sed 's/color=//'| sed 's/label="//'| sed 's/"//g'| sed 's/]//'| sed 's/\\n/__/' |sed 's/^\t//'> T-PRL-S-binarized_nodes.txt


nodes = sys.argv[1]
edges = sys.argv[2]
coldict = {}
nodedict = {}

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
        if nodeid not in nodedict:
            nodedict[nodeid] = name
        # ignore redundant labels of false polytomy resolutions
        if name.count('^') < 2 or name.split("^")[2]=="0":
            #print("label",name,tissue)
            print("label",col,name)
            

with open(edges,'r') as e:
    for l in e:
        l = l.strip()
        l = l.split(" ")
        parent = l[0]
        child = l[1]
        parent_name = nodedict[parent]
        child_name = nodedict[child]
        # ignore redundant parent-child relationships caused by false polytomy resolutions
        if not parent_name.count('^') > 1 and not child_name.count('^') > 1:
            #if parent_name.count('^') > 1:
            #    parent_name = parent_name.split("^")[0] + "^" + parent_name.split("^")[1] + "^0"
            print("tree",parent_name,child_name)
        elif parent_name.count('^') > 1 and not child_name.count('^') > 1:
            if not parent_name.split("^")[2]=="0":
                parent_name = parent_name.split("^")[0] + "^" + parent_name.split("^")[1] + "^0"
            print("tree",parent_name,child_name)
        elif not parent_name.count('^') > 1 and child_name.count('^') > 1:
            if child_name.split("^")[2]=="0":
                print("tree",parent_name,child_name)
        #elif parent_name.count('^') > 1 and child_name.count('^') > 1:
            #print("ignore case:",parent_name,child_name)


# nodes
#277 3 ASV030_PRL\nPRL
#276 1 ASV222_LVM\nLVM

# parent child
#274 276
#274 275
