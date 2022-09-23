import sys
from collections import defaultdict
from collections import Counter

asv_dict = defaultdict(lambda: defaultdict(set))
asv_per_dict = defaultdict(list)
asv_done = []
with open(sys.argv[1],'r') as infile:
    header = next(infile)
    for l in infile:
        l = l.strip()
        l = l.split(",")
        asv = l[0]
        mutation = l[25]
        cp = l[29]
        asv_count = int(l[3])
        asv_per = float(l[5])
        asv_dict[cp][asv].add(mutation)
        if asv not in asv_done:
            asv_per_dict[cp].append(asv_per)
            asv_done.append(asv)

print(asv_dict)
for k,v in asv_dict.items():
    asv_list = []
    for amplicon,mutations in v.items():
        asv_list.append(sorted(list(mutations)))
    mutation_counts = Counter([tuple(i) for i in asv_list])
    per_uniq = len(mutation_counts) / len(v) 
    print(k,per_uniq,len(mutation_counts),len(v),sorted(asv_per_dict[k]))
