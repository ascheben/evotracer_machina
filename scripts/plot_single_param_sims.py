import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

### This script takes in simulation output where all variables were 
### consant except for one variable in which a range was tested
### The first command line argument is the filepath to csv output from parallel_sim.sh
### The second command line argument is the column header name for the variable of interest in the csv file
### The third command line argument is the filepath to same output directory from parallel_sim.sh 
  
def plot_csv(csv_file, x_column, outdir):
    data = pd.read_csv(csv_file)

    ### Can further stratify the data based on other parameters as constants if needed.
    #data = data[data['mutrate']==0.05]
    #data = data[data['num_samples']==100]
    data_rare = data[data['migration_matrix']=='data/rare_migration_prob_matrix.csv']
    data_true = data[data['migration_matrix']=='data/true_migration_prob_matrix.csv']
    data_equal = data[data['migration_matrix']=='data/equal_migration_prob_matrix.csv']

    data['num_sites']=(np.round(np.array(data[x_column]/data['mut_per_site']), 0))

    #data_rare = data_rare.sort_values(x_column)
    #data_true = data_true.sort_values(x_column)
    #data_equal = data_equal.sort_values(x_column)
    ci=95
    plt.figure(figsize=(12, 8))
    #sns.lineplot(x=x_column, y='proportion', data=data_rare, errorbar=('ci', ci), err_style='band', marker='o', color='red', label='Rare migration matrix')
    #sns.lineplot(x=x_column, y='proportion', data=data_true, errorbar=('ci', ci), err_style='band', marker='o', color='blue', label='True migration matrix')
    #sns.lineplot(x=x_column, y='proportion', data=data_equal, errorbar=('ci', ci), err_style='band', marker='o', color='green', label='Equal migration matrix')
    #sns.barplot(x=x_column, y='proportion', data=data, hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'], dodge=True, errorbar=('ci', ci), palette=['red', 'blue', 'green'], capsize=0.1)
    #sns.barplot(x='num_sites', y='proportion', data=data, hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'], dodge=True, errorbar=('ci', ci), palette=['red', 'blue', 'green'], capsize=0.1)
    #sns.scatterplot(x=x_column, y='proportion', data=data_rare, marker='o', color='red', label='Rare migration matrix')
    #sns.scatterplot(x=x_column, y='proportion', data=data_true, marker='o', color='blue', label='True migration matrix')
    #sns.scatterplot(x=x_column, y='proportion', data=data_equal, marker='o', color='green', label='Equal migration matrix')
    #sns.boxplot(x=x_column, y='proportion', data=data_rare, color='red', width=0.2, showfliers=False, dodge=True)
    #sns.boxplot(x=x_column, y='proportion', data=data_true, color='blue', width=0.2, showfliers=False, dodge=True)
    #sns.boxplot(x=x_column, y='proportion', data=data_equal, color='green', width=0.2, showfliers=False, dodge=True)
    #sns.boxplot(x=x_column, y='proportion', hue='migration_matrix', data=data, width=0.5, palette=['red', 'blue', 'green'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'], showfliers=False, dodge=True)
    #sns.swarmplot(dodge=True, size=1, x=x_column, y='proportion', hue='migration_matrix', data=data, palette=['red', 'blue', 'green'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    #sns.stripplot(dodge=True, size=2, x=x_column, y='proportion', hue='migration_matrix', data=data, palette=['red', 'blue', 'green'], hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    sns.barplot(x='num_sites', y='proportion', data=data, dodge=True, errorbar=('ci', ci), capsize=0.1, palette=['red', 'blue', 'green'], hue='migration_matrix', hue_order=['data/rare_migration_prob_matrix.csv', 'data/true_migration_prob_matrix.csv', 'data/equal_migration_prob_matrix.csv'])
    plt.xlabel('num_sites')
    #plt.xticks(rotation=-80)
    #plt.ylim(0, 1)
    plt.ylabel('Proportion of true migrations inferred')
    plt.title('95 CI')
    #plt.legend(['Rare migration matrix', 'True migration matrix', 'Equal migration matrix'], loc='upper left', bbox_to_anchor=(0, 1))
    #plt.legend(loc='upper left', bbox_to_anchor=(0, 1))
    plt.tight_layout()
    plt.savefig(f'{outdir}numSites_vs_proportion_barplot.png')
    plt.show()

csv=sys.argv[1]
single_param_name=sys.argv[2]
out_dir=sys.argv[3]

plot_csv(csv, single_param_name, out_dir)
