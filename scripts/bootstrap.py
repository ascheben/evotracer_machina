import pandas as pd
import numpy as np
import cassiopeia as cas

df = pd.read_csv('../data/MMUS1469_matrix.csv',index_col=0)

# exclude trunkal mutations from bootstrapping and add them back later
#trunkal = ["i_147_1nts","i_17_1nts","d_15_37nts","d_19_10nts","d_18_2nts","d_17_1nts","d_13_8nts","d_19_1nts","d_8_26nts","d_18_9nts","d_11_10nts","d_18_3nts","d_147_1nts","d_18_25nts","d_139_11nts","d_18_4nts","i_17_2nts","d_146_2nts","d_148_25nts","d_20_26nts","i_147_2nts","d_6_13nts","d_18_8nts","d_147_30nts","d_128_52nts","d_13_4nts","d_147_4nts","d_18_28nts","d_153_55nts","d_129_22nts","d_16_3nts","d_12_6nts","d_147_10nts","d_149_1nts","d_148_3nts","d_148_10nts","d_177_52nts","d_17_5nts","d_146_3nts","d_13_10nts","d_147_5nts","d_17_27nts","d_14_4nts","i_14_5nts","d_8_28nts","i_149_1nts","d_9_40nts","d_45_3nts","d_17_11nts","d_140_26nts","d_136_39nts","d_124_26nts","d_10_54nts","d_141_12nts","d_147_28nts","d_137_20nts","d_144_49nts","d_18_24nts","i_43_1nts","d_16_2nts","d_145_2nts","d_149_10nts","d_166_10nts","i_20_5nts"]
#trunkal_cols = df[df.columns & trunkal]
#df = df.drop(trunkal, axis=1)
#cas_tree = cas.data.CassiopeiaTree(character_matrix=df)
cas_boot = cas.data.sample_bootstrap_character_matrices(character_matrix=df,num_bootstraps=100)
vanilla_greedy = cas.solver.VanillaGreedySolver()

for cas_mat in cas_boot:
    # add the trunkal mutations back to the matrix
    #new_mat = cas_mat[0].join(trunkal_cols)
    #cas_tree = cas.data.CassiopeiaTree(character_matrix=new_mat)
    cas_tree = cas.data.CassiopeiaTree(character_matrix=df)
    vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)
    cas_tree.get_newick(record_branch_lengths = False)
    cas_tree.collapse_mutationless_edges(infer_ancestral_characters=True)
    ntree = cas_tree.get_tree_topology()
    for n, nbrs in ntree.adjacency():
        if len(nbrs) > 0:
            subdict = nbrs
            nmut_parent = sum(cas_tree.get_character_states(n))
            for child in subdict:
                nmut_child = sum(cas_tree.get_character_states(child))
                cas_tree.set_branch_length(n,child,float((nmut_child - nmut_parent)))
    
    print(cas_tree.get_newick(record_branch_lengths = True))
