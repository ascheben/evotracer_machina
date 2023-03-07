import sys
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree

tabular_tree = []
with open(sys.argv[1],'r') as tabtree:
    for l in tabtree:
        l = l.strip().split(" ")
        tabular_tree.append((l[0],l[1],1))
tree = Tree.from_parent_child_table(tabular_tree)
outnwk = sys.argv[1].split(".")[0] + ".nwk"
print(tree.write(format=1,outfile=outnwk))
outpdf = sys.argv[1].split(".")[0] + ".pdf"
tree.render(outpdf)
