import sys
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree

t = Tree(sys.argv[1])
outpdf = sys.argv[1].split(".")[0] + ".pdf"

t.render(outpdf)

