import sys
import numpy as np 

ref_cut_sites = [17, 43, 69, 95, 121, 147, 173, 199, 225, 251]
ref_border_sites = [1, 26, 52, 78, 104, 130, 156, 182, 208, 234]
def closest(lst, K):
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx]

with open(sys.argv[1],'r') as asv:
    next(asv)
    for l in asv:
        l = l.strip().split('\t')
        seqname = l[0]
        #indeltype = l[1]
        tag = l[7]
        # nearest cut site to leftmost position of indel, favoring within-indel sites
        # is there one cut site within indel?
        indel_type = tag.split('_')[0]
        indel_start = int(tag.split('_')[1])
        indel_len = int(tag.split('_')[2].replace('nts',''))
        indel_end = indel_start + indel_len - 1
        # more than 1 cut site?
        overlap_cut_sites = []
        for site in ref_cut_sites:
            if site >= indel_start and site<= indel_end:
                overlap_cut_sites.append(site)
        # is cut site ambiguous?
        if len(overlap_cut_sites) > 1:
            dist = "NA"
            assign_cut_site = "NA"
        elif len(overlap_cut_sites) == 1:
            dist = 0
            assign_cut_site = overlap_cut_sites[0]
        else:
            nearest = closest(ref_cut_sites,indel_start)
            dist = nearest - indel_start
            assign_cut_site = nearest
        print(*[seqname,indel_type,indel_len,indel_start,indel_end,dist,assign_cut_site],sep='\t')
        # this is the case if the indel overlaps multiple cut sites
        # no within indel sites
        # what is closest cut site to leftmost pos


#SEQ11072        d       144     150     34      TCAACC      ------    d_144_6nts     141     153     147     2
#SEQ11071        d       144     177     34      TCAACCAGGGTCAGATACACCCACGTTCAACCGT      ----------------------------------    d_144_34nts     141     180     c(147, 173)     2
