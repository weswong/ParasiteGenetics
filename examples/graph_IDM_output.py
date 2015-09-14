# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import seaborn as sns
import json
import sys
import os
import numpy as np
import glob
from scipy import stats
import scipy as sp

# <codecell>

def extract_data(file):
    data= json.load(open(file))
    time = np.asarray(data['tsteps'])
    f_polygenomic = np.asarray(data['f_polygenomic']['Population #1'])
    f_infected = np.asarray(data['f_infected']['Population #1'])
    return time, f_polygenomic, f_infected

# <codecell>

def run(cwd):
    data_dict = {'time':[], 'f_polygenomic':[], 'f_infected':[]}
    for dir in glob.glob('simulations_*'):
        working_dir =  dir
        os.chdir(working_dir)
        time, f_polygenomic, f_infected = extract_data('PopulationInfectionReport.json')
        data_dict['time'].append( time )
        data_dict['f_polygenomic'].append( f_polygenomic )
        data_dict['f_infected'].append( f_infected )
        os.chdir(cwd)
    return data_dict

    

# <codecell>

class graph_stats():
    def __init__(self, data_array):
        self.mean = np.mean(data_array, axis=0)
        self.std = np.std(data_array, axis=0)
        
        self.confidence= [stats.norm.interval(0.95, loc=avge, scale=self.std[i]) for i,avge in enumerate(self.mean)]

# <codecell>

def make_graph(data_dict):
    
    f_infected = graph_stats(data_dict['f_infected'])
    f_polygenomic = graph_stats(data_dict['f_polygenomic'])
    
    #f_infected
    plt.subplot(2,1,1)
    plt.title('Proportion of people infected and proportion of mixed infections')
    plot(time, f_infected.mean, label='f_infected', color = 'green')
    plot(time, [x[0] for x in f_infected.confidence], color = 'green', alpha=0.5)
    plot(time, [x[1] for x in f_infected.confidence], color = 'green', alpha=0.5)
    plt.fill_between(time, [x[0] for x in f_infected.confidence], [x[1] for x in f_infected.confidence], color='green', alpha='0.3')
    #real data
    #plt.scatter([365* (x + 5) for x in range(0,8)], ['0.114386', '0.0566608', '0.00882555', '0.00712973', 
    #                                                 '0.0278191', '0.0376105', '0.0216132', '0.0533142'],
    #            color = 'green')
    plt.legend(loc='lower left', fontsize=15)
    plt.subplot(2,1,2)
    plot(time, f_polygenomic.mean, label='f_polygenomic', color = 'blue')
    plot(time, [x[0] for x in f_polygenomic.confidence], color = 'blue', alpha =0.5)
    plot(time, [x[1] for x in f_polygenomic.confidence], color = 'blue', alpha =0.5)
    plt.fill_between(time, [x[0] for x in f_polygenomic.confidence], [x[1] for x in f_polygenomic.confidence], color='blue', alpha='0.3')
    #real data
    #plt.scatter([365* (x + 5) for x in range(0,8)], ['0.312977', '0.334951', '0.119266', '0.177083', 
    #                                                 '0.206349', '0.188406', '0.410448', '0.251592'],
    #            color = 'blue')
    
    
    plt.legend(loc='lower left', fontsize=15)
    
    plt.savefig('graph.png')

# <codecell>

if __name__ == '__main__':
    sns.set_style(style='white', rc=None)
    cwd = os.getcwd()
    run(cwd)
    make_graph(data_dict)

# <codecell>


