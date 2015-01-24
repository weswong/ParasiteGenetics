import itertools

import numpy as np

from parameter_sweep import parameter_sweep
from visualize_sweep import visualize

'''
"R0_initial" : 2.7
"R0_reduction_year" : 7,
"R0_final" :  1.95,
"R0_plateau_year" : 10,
"R0_rebound_year" : 12,
"R0_rebound" : 2.35,
"rates" : { "expire" :  0.022,
            "import" :  0.01 },
"seasonality" : { "minimum" : 0.18,
                "coefficient" : 2 },
"n_populations" : 1,
"initial_parasites" : 1200,
"initial_unique" : 1000,
"n_humans" : 1350
'''

param_ranges = {
  #'R0_rebound_year' : np.arange(10.5,14.5,0.5),
  #'R0_rebound'      : np.arange(1.95,2.85,0.1),
  'R0_plateau_year' : np.arange(8.0,12.0,0.5),
  'R0_final'        : np.arange(1.55,2.35,0.1),
}

param_names  = param_ranges.keys()
param_combos = itertools.product(*param_ranges.values())
work_dir     = 'simulations/'+'_'.join(param_names)

parameter_sweep(param_names,param_combos,work_dir,nseeds=5,num_cpus=3)

visualize(work_dir,param_names)