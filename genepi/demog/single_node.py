import logging
log = logging.getLogger(__name__)

import demog

def init():

    populations = {
        'Population #1' : {
            'n_humans' : 1350,
            'n_infections' : 800,
            'coi' : demog.init_coi_discrete_distribution(n=5, alpha=0.5, beta=8),
            'r0_params': {'r0' : [1, 1, 1, 1.5, 1, 1],'day' : [0, 5*365,10*365,15*365,20*365,30*365]},
            'annual_cycle': demog.annual_cycle(0.15,1e-2,8),
                            
        }
    }

    return populations
