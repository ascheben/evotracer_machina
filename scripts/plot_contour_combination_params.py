
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate


def plot_contour(data_filepath, var1, var2, result_var, output_dir):
    # Load data from CSV file
    raw_data = pd.read_csv(f'{data_filepath}')
    data = raw_data[raw_data['migration_matrix']=='data/equal_migration_prob_matrix.csv']    ### change this if stratifying on another categorical parameter

    # Extract variables from dataframe
    x = data[f'{var1}'].values
    y = data[f'{var2}'].values
    z = data[f'{result_var}'].values

    # Define grid for contour plot
    xgrid = np.linspace(min(x), max(x), 100)
    ygrid = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xgrid, ygrid)

    # Interpolate effect variable onto grid
    Z = scipy.interpolate.griddata((x, y), z, (X, Y), method='linear')

    # Define color scale and colorbar legend
    levels = np.linspace(min(z), max(z), 10)
    #levels = np.linspace(0, 1) ### use this line to set 0 and 1 limits to normalize the proportion values
    cmap = plt.get_cmap('coolwarm')
    norm = plt.Normalize(vmin=min(levels), vmax=max(levels))
    fig, ax = plt.subplots()
    CS = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, norm=norm)
    cbar = fig.colorbar(CS)
    cbar.ax.set_ylabel(f'{result_var}')
    
    CS_lines = ax.contour(X, Y, Z, colors='k', levels=levels)
    #ax.clabel(CS_lines, inline=1, fontsize=10)
    ax.set_title('Equal migration matrix')   ### Fix for categorical parameter name
    ax.set_xlabel(f'{var1}')
    ax.set_ylabel(f'{var2}')
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{output_dir}/contour_equalMM2_{var1}_{var2}_{result_var}.png')     ### Can specify the name further for categorical stratification


### Takes in filepath to data csv, x variable, y variable, the z measured variable, and the desired output directory for saving results.
data_path = sys.argv[1]
ind_variable1 = sys.argv[2]
ind_variable2 = sys.argv[3]
measured_variable = sys.argv[4] 
outdir = sys.argv[5]

plot_contour(data_path, ind_variable1, ind_variable2, measured_variable, outdir)

