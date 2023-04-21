import sys
import pandas as pd
import matplotlib.pyplot as plt

### This script takes in simulation output where all variables were 
### consant except for one variable in which a range was tested
### The first command line argument is the filepath to csv output from parallel_sim.sh
### The second command line argument is the column header name for the variable of interest in the csv file
### The third command line argument is the filepath to same output directory from parallel_sim.sh 
  
def plot_csv(csv_file, x_column, outdir):
    data = pd.read_csv(csv_file)
    data = data.sort_values(x_column)
    x = data[x_column]
    y = data['proportion']
    plt.plot(x, y, marker='o')
    plt.xlabel(x_column)
    #plt.xticks(rotation=20,fontsize=10)
    plt.ylim(0, 1)
    plt.ylabel('Proportion of true migrations inferred')
    plt.title('Proportion vs ' + x_column)
    plt.tight_layout()
    plt.savefig(f'{outdir}proportion_vs_{x_column}.png')
    plt.show()

csv=sys.argv[1]
single_param_name=sys.argv[2]
out_dir=sys.argv[3]

plot_csv(csv, single_param_name, out_dir)
