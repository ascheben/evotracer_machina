import pandas as pd
import numpy as np
import sys

def calculate_f1_score(sim_mut, result_mut):
    common_result_mut = np.intersect1d(sim_mut, result_mut)

    true_positives = len(common_result_mut)
    false_positives = len(result_mut) - true_positives
    false_negatives = len(sim_mut) - true_positives

    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    f1_score = 2 * (precision * recall) / (precision + recall)

    return precision, recall, f1_score

def calculate_f1_score_flexible(sim_mut, result_mut):
    # Find the perfect matches between sim_mut and result_mut
    common_result_mut = np.intersect1d(sim_mut, result_mut)
    # Create an array of unmatched mutations in result_mut and sim_mut
    result_unmatched = np.setdiff1d(result_mut, common_result_mut)
    sim_unmatched = np.setdiff1d(sim_mut, common_result_mut)
    # Filter out deletions and mutations shorter than or equal to 4bp
    filtered_result_unmatched = [m for m in result_unmatched if '-' not in m and len(m) > 3]
    filtered_sim_unmatched = [m for m in sim_unmatched if '-' not in m and len(m) > 3]

    true_positives = len(common_result_mut)

    #add to true positives if match is only 1bp off on the end for insertions
    for sim_mutation in filtered_sim_unmatched:
        match_found = False

        for result_mutation in filtered_result_unmatched:
            # Check if the mutation matches with 1bp removed from the left
            shifted_result_mutation = result_mutation[1:]
            shifted_sim_mutation = sim_mutation[:-1]
            if np.array_equal(shifted_sim_mutation, shifted_result_mutation):
                true_positives += 1
                match_found = True
                continue
            shifted_result_mutation = result_mutation[:-1]
            shifted_sim_mutation = sim_mutation[1:]
            if np.array_equal(shifted_sim_mutation, shifted_result_mutation):
                true_positives += 1
                match_found = True
                continue

    false_positives = len(result_mut) - true_positives
    false_negatives = len(sim_mut) - true_positives

    # Calculate precision, recall, and F1 score
    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    f1_score = 2 * (precision * recall) / (precision + recall)

    return precision, recall, f1_score


sim_name= sys.argv[1]
mutrate=sys.argv[2]
n_samp=sys.argv[3]
sim_tsv=sys.argv[4]
amp_csv = sys.argv[5]
output_dir=sys.argv[6]
amp_read_count_csv = sys.argv[7]
evo_csv= sys.argv[8]
cresso_txt = sys.argv[9]

#sim_tsv="sim_results_evotracerAndAmplicanSampleData/out_simulator_simmid/simmid_mutations.tsv"
#amp_csv="sim_results_evotracerAndAmplicanSampleData/out_amplican_simmid/alignments/events_filtered_shifted.csv"
#evo_csv="sim_results_evotracerAndAmplicanSampleData/out_evotracer_parallelize_simmid/phylogeny_analysis/phylogeny_del_ins/asv_stat.csv"

sim_df = pd.read_csv(sim_tsv, sep='\t', header=None)
if amp_csv != "NA":
    amp_df = pd.read_csv(amp_csv)
    amp_reads_df = pd.read_csv(amp_read_count_csv)
if evo_csv != "NA":
    evo_df = pd.read_csv(evo_csv)
if cresso_txt != "NA":
    cresso_df = pd.read_csv(cresso_txt, delim_whitespace=True)

## obtain sim info
sim_mut = sim_df.iloc[:,1].unique()
sim_mut_string = ' '.join(sim_mut.astype(str))
sim_mut_count = len(sim_mut)

## obtain EvoTraceR info
if evo_csv != "NA":
    evo_df = evo_df.iloc[:-3,:] ##remove last three rows corresponding to unmutated barcode for 3 tissues by default
    evo_mut = evo_df['alt_seq'].unique()
    evo_mut_string = ' '.join(evo_mut.astype(str))
    evo_mut_count = len(evo_mut)
    evo_precision, evo_recall, evo_f1 = calculate_f1_score_flexible(sim_mut, evo_mut)
else:
    evo_mut_string = "NA"
    evo_mut_count = "NA"
    evo_precision = "NA"
    evo_recall = "NA"
    evo_f1 = "NA"

## obtain ampliCan info
if amp_csv != "NA":
    total_reads = amp_reads_df['filtered_read_count']
    amp_df.loc[amp_df['type'] == 'deletion', 'replacement'] = amp_df.loc[amp_df['type'] == 'deletion', 'width'].apply(lambda x: '-' * x) ##replace "" in replacement with dashes the length of the width column for deletions to match evo format
    amp_df = amp_df.loc[amp_df['type'].isin(['insertion','deletion'])]
    collapsed_df = amp_df.groupby('replacement').agg({'counts': 'sum'}).reset_index()
    #amp_mut = amp_df['replacement'].unique()
    collapsed_df['freq'] = collapsed_df['counts'].apply(lambda x: x/total_reads)
    collapsed_df = collapsed_df.loc[collapsed_df['freq'] > 0.01]  ## default frequency for amplicanPipeline to define sequencing errors
    amp_mut = collapsed_df['replacement']
    amp_mut_string = ' '.join(amp_mut.astype(str))
    amp_mut_count = len(amp_mut)
    amp_precision, amp_recall, amp_f1 = calculate_f1_score_flexible(sim_mut, amp_mut)
else:
    amp_mut_string = "NA"
    amp_mut_count = "NA"
    amp_precision = "NA"
    amp_recall = "NA"
    amp_f1 = "NA"

## obtain CRISPResso2 info
if cresso_txt != "NA":
    ref_seq = cresso_df.loc[0,'Aligned_Sequence']
    cresso_df = cresso_df.loc[cresso_df['Unedited'] == False]
    # call mutations
    for index, row in cresso_df.iterrows():
        size = 0
        if row['n_deleted'] != 0:
            size = row['n_deleted']
        if row['n_inserted'] != 0:
            size = row['n_inserted']
        end = 20 + size
        cresso_df.at[index,'mutation'] = row['Aligned_Sequence'][20:end]
    cresso_df = cresso_df.groupby('mutation', as_index=False).sum()
    cresso_df = cresso_df.sort_values('%Reads', ascending=False)
    cresso_df = cresso_df.loc[cresso_df['%Reads'] > 0.2]  ## default freq cutoff for CRISPRessso2 documentation
    cresso_mut = cresso_df['mutation']
    cresso_mut_string = ' '.join(cresso_mut.astype(str))
    cresso_mut_count = len(cresso_mut)
    cresso_precision, cresso_recall, cresso_f1 = calculate_f1_score_flexible(sim_mut, cresso_mut)
else:
    cresso_mut_string = "NA"
    cresso_mut_count = "NA"
    cresso_precision = "NA"
    cresso_recall = "NA"
    cresso_f1 = "NA"
    

# make result dataframe
# if evo_csv != "NA":
cols = ['sim_name','mutrate','num_samples','sim_mutations','count_sim_mutations','evotracer_mutations','count_evotracer_mutations','evo_precision','evo_recall','evo_f1','amplican_mutations','count_amplican_mutations','amp_precision','amp_recall','amp_f1','crispresso2_mutations','count_crispresso2_mutations','crispresso2_precision','crispresso2_recall','crispresso2_f1']
d = [sim_name,mutrate,n_samp,sim_mut_string,sim_mut_count,evo_mut_string,evo_mut_count,evo_precision,evo_recall,evo_f1,amp_mut_string,amp_mut_count,amp_precision,amp_recall,amp_f1,cresso_mut_string,cresso_mut_count,cresso_precision,cresso_recall,cresso_f1]

# else:
#     cols = ['sim_name','mutrate','num_samples','sim_mutations','count_sim_mutations','amplican_mutations','count_amplican_mutations','amp_precision']
#     d = [sim_name,mutrate,n_samp,sim_mut_string,sim_mut_count,amp_mut_string,amp_mut_count,precision_amp_mut]

result_df = pd.DataFrame(
                        data=[d],
                        columns=cols
)

output_file = f'{output_dir}compare_amplican_evotracer_crispresso_{sim_name}.csv'
result_df.to_csv(output_file, index=False)

