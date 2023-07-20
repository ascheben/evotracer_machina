import sys
import re

#grep label G-PRL-S-binarized.dot| tr ',' ' '| cut -d' ' -f1,6| sed 's/\t//g'| sed 's/color=//'| sed 's/label=//'| sed 's/"//g'| sed 's/]//' | sed 's/^\t//'> G-PRL-S-binarized_labels.txt
#grep "\->" G-PRL-S-binarized.dot| cut -d' ' -f1,3|sed 's/^\t//' > G-PRL-S-binarized_parent_child.txt

tissues = sys.argv[1]
migrations = sys.argv[2]
tisdict = {}

with open(tissues,'r') as n:
    for l in n:
        l = l.strip()
        l = l.split(" ")
        t_id = l[0]
        tissue = l[1]
        if t_id not in tisdict:
            tisdict[t_id] = tissue

with open(migrations,'r') as e:
    for l in e:
        l = l.strip()
        l = l.split(" ")
        parent = l[0]
        child = l[1]
        parent_name = tisdict[parent]
        child_name = tisdict[child]
        print(parent_name,child_name)


# colors
#2 PRL
#3 RBL
#0 LVM
#1 LVR

#migrations
#2 1
#2 1
