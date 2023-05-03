import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

### This script takes in simulation output where all variables were 
### consant except for one variable in which a range was tested
### The first command line argument is the filepath to csv output from parallel_sim.sh
### The second command line argument is the column header name for the variable of interest in the csv file
### The third command line argument is the filepath to same output directory from parallel_sim.sh 
  
def plot_csv(csv_file, x_column, outdir):
    data = pd.read_csv(csv_file)
    #data=data.fillna(0)

    ### Data import to be used for scipy result plotting
    #data = pd.read_csv('SIM_RESULTS_STEPHEN/7_constantMutrate_scipy_results/data_constantMutrate_scipy_results.csv',usecols=range(1,102,1))
    #data_melt = pd.melt(data, id_vars='avg_mutrate')

    ### Can further stratify the data based on other parameters as constants if needed.
    #data = data[data['mutrate']==0.05]
    #data = data[data['num_samples']==100]
    #data_rare = data[data['migration_matrix']=='data/rare_migration_prob_matrix.csv']
    #data_true = data[data['migration_matrix']=='data/true_migration_prob_matrix.csv']
    #data_moderate = data[data['migration_matrix']=='data/moderate_migration_prob_matrix.csv']
    #data_high = data[data['migration_matrix']=='data/high_migration_prob_matrix.csv']
    #data_equal = data[data['migration_matrix']=='data/equal_migration_prob_matrix.csv']

    #data['average_mutrate'] = data['mutrate'].apply(lambda x: sum(map(float, x.split()))/len(x.split()))

    #def avg_greater_than_zero(x):
    #    values = list(map(float, x.split()))
    #    values = [v for v in values if v > 0]
    #    return sum(values)/len(values) if len(values) > 0 else 0

    #def avg_greater_than_zero(x):
     #   values = list(map(float, x.split()))
      #  if all(v > 0 for v in values):
       #     return sum(values)/len(values)
        #else:
         #   return 0

    #data['average_mutrate'] = data['mutrate'].apply(avg_greater_than_zero)
    #data = data[data['average_mutrate'] != 0]


    #def all_ten_active(x):
    #    values = list(map(float, x.split()))
    #    values = [v for v in values if v > 0]
    #    return True if len(values) == 10 else False

    #data['all_ten'] = data['mutrate'].apply(all_ten_active)
    #data = data[data['all_ten']==True]
    
    #mutrates=data['mutrate'].unique()
    #data = data[data['mutrate'].isin(mutrates[0:10])]

    data[x_column]=(np.round(np.array(data[x_column]), 1))

    #data['average_mutrate']=(np.round(np.array(data['average_mutrate']), 2))

    #data['avg_proportion_mut_sites']=(np.round(np.array(data[x_column]), 0))
    #data[x_column] = data[x_column].round(1)
    #data_saturation=data[data[x_column].astype('float64') > 46]

    #data_melt = pd.melt(data, id_vars=[x_column], value_vars=['true_migrations', 'inferred_migrations'], 
    #                    var_name='migration_type', value_name='migrations')

    #data_rare = data_rare.sort_values(x_column)
    #data_true = data_true.sort_values(x_column)
    #data_equal = data_equal.sort_values(x_column)
    ci=95
    plt.figure(figsize=(12, 8))
    #sns.lineplot(x=x_column, y='proportion', data=data_rare, errorbar=('ci', ci), err_style='band', marker='o', color='red', label='Rare migration matrix')
    #sns.lineplot(x=x_column, y='proportion', data=data_true, errorbar=('ci', ci), err_style='band', marker='o', color='blue', label='True migration matrix')
    #sns.lineplot(x=x_column, y='proportion', data=data_equal, errorbar=('ci', ci), err_style='band', marker='o', color='green', label='Equal migration matrix')
    sns.lineplot(x=x_column, y='num_uniq_mutations', data=data, errorbar=('ci', ci), err_style='band', marker='o')
                #,palette=['red', 'blue', 'green', 'black', 'orange'], hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/moderate_migration_prob_matrix.csv', 'data/high_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    #sns.scatterplot(x=x_column, y='proportion', data=data, size=2, palette=['red', 'blue', 'green', 'black', 'orange'], hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/moderate_migration_prob_matrix.csv', 'data/high_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    #sns.barplot(x=x_column, y='proportion', data=data, hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'], dodge=True, errorbar=('ci', ci), palette=['red', 'blue', 'green'], capsize=0.1)
    #sns.barplot(x='num_sites', y='proportion', data=data, hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'], dodge=True, errorbar=('ci', ci), palette=['red', 'blue', 'green'], capsize=0.1)
    #sns.scatterplot(x=x_column, y='proportion', data=data_rare, marker='o', color='red', label='Rare migration matrix')
    #sns.scatterplot(x=x_column, y='proportion', data=data_true, marker='o', color='blue', label='True migration matrix')
    #sns.scatterplot(x=x_column, y='proportion', data=data_equal, marker='o', color='green', label='Equal migration matrix')
    #sns.boxplot(x=x_column, y='proportion', data=data_rare, color='red', width=0.2, showfliers=False, dodge=True)
    #sns.boxplot(x=x_column, y='proportion', data=data_true, color='blue', width=0.2, showfliers=False, dodge=True)
    #sns.boxplot(x=x_column, y='proportion', data=data_equal, color='green', width=0.2, showfliers=False, dodge=True)
    #sns.boxplot(x=x_column, y='proportion', hue='migration_matrix', data=data, width=1, palette=['red', 'blue', 'green', 'black', 'orange'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/moderate_migration_prob_matrix.csv', 'data/high_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'], showfliers=False, dodge=True)
    #sns.swarmplot(dodge=True, size=1, x=x_column, y='proportion', hue='migration_matrix', data=data, palette=['red', 'blue', 'green', 'grey'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    #sns.stripplot(dodge=True, size=2, x=x_column, y='proportion', hue='migration_matrix', data=data, palette=['red', 'blue', 'green', 'grey'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    #sns.barplot(x=x_column, y='proportion', data=data, dodge=True, errorbar=('ci', ci), capsize=0.1, palette=['red', 'blue', 'green', 'grey'], hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/moderate_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    

    #sns.boxplot(x=x_column, y='migrations', data=data_melt, hue='migration_type', linewidth=1.5, width=0.4, palette=['black', 'red'], showfliers=False, dodge=True)
    #sns.stripplot(x=x_column, y='migrations', data=data_melt, hue='migration_type', size=2, palette=['red', 'blue'], dodge=True, jitter=True)

    #sns.stripplot(x='mutrate', y='proportion', hue='migration_matrix', data=data_saturation, jitter=True, palette=['red', 'blue', 'green', 'orange'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/moderate_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'],dodge=True)

    # sns.regplot(x='average_mutrate', y='proportion', data=data, color='blue')
    # statistic,pvalue = stats.pearsonr(x=data['average_mutrate'], y=data['proportion'])
    # plt.annotate(f"Pearson Corr: {statistic:.2f}\nP value: {pvalue}", 
    #          xy=(0.6, 0.9), 
    #          xycoords='axes fraction',
    #          fontsize=12)

    #sns.scatterplot(x='average_mutrate', y='avg_proportion_mut_sites', data=data, s=20, color='k', legend=False)
    #sns.lineplot(x='avg_mutrate', y='value', data=data_melt, color='k', errorbar=('ci', ci), err_style='band', marker='o')
    #sns.boxplot(x=x_column, y='avg_mutation_age', data=data, width=0.6, color='grey', showfliers=False, fliersize=3, dodge=True)
    
    ### Used below to annotate scipy lineplot
    #plt.annotate(f"Scipy.optimize solution:\navg_mutrate = 0.126398\nproportion = 0.60813", 
    #            xy=(0.68, 0.85), 
    #            xycoords='axes fraction',
    #            fontsize=18,
    #            bbox = dict(boxstyle="round,pad=0.3", fc="none", ec="r", lw=2))
    #plt.scatter(0.12639846516042785, 0.6081340000000001, s=3000, facecolors='none', linewidth=2, edgecolors='r')



    plt.xlabel(x_column)
    #plt.xticks(ticks=range(0,1,0.1), labels=range(0,1,0.1))
    #plt.xticks(rotation=-40, ha='left')
    #plt.xlim(0.12635,0.12643)
    plt.ylabel('num_uniq_mutations')
    plt.title('')
    #plt.legend(['Rare migration matrix', 'True migration matrix', 'Equal migration matrix'], loc='upper left', bbox_to_anchor=(0, 1))
    plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
    plt.tight_layout()
    plt.savefig(f'{out_dir}{x_column}_vs_numUniqMutations_lineplot.png')
    #plt.show()

csv=sys.argv[1]
single_param_name=sys.argv[2]
out_dir=sys.argv[3]

plot_csv(csv, single_param_name, out_dir)
