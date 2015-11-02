import logging
log = logging.getLogger(__name__)

import demog

def init():

    populations = {
        'Population #1' : {
            'n_humans' : 16000,
            'n_infections' : 1600,
            'coi' : demog.init_coi_discrete_distribution(n=5, alpha=0.5, beta=8),
            'r0_params': {'r0' : [1, 1, 0.30, 1.2, 1.0],
                          'day' : [0,5*365,7*365, 8*365, 12*365]},
            'annual_cycle': lambda t: 0.0365
            #demog.annual_cycle(0.15,1e-2,8),
                            
        }
    }

    return populations
