import sys

#ASV003,d_15_37nts,CP04
#asv_stat_names_mut_cp.csv
cp_dict = {}
indel_dict = {}

with open(sys.argv[1],'r') as asv:
    for l in asv:
        l = l.strip()
        l = l.split(",")
        name = l[0]
        indel = l[1]
        cp = l[2]
        if cp not in cp_dict:
            cp_dict[cp] = set()
            cp_dict[cp].add(name)
        else:
            cp_dict[cp].add(name)
        if name not in indel_dict:
            indel_dict[name] = set()
            indel_dict[name].add(indel)
        else:
            indel_dict[name].add(indel)
for c,name_set in cp_dict.items():
    # find all indels in cp
    cp_indels = set()
    for name in name_set:
        cp_indels.update(list(indel_dict[name]))
    for indel in cp_indels:
        all_have = True
        for name in name_set:
            if indel not in indel_dict[name]:
                all_have = False
        if all_have:
            #print(f'Truncal mutation of {c} is {indel}')
            print(c,indel)
