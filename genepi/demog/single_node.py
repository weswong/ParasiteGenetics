import logging
log = logging.getLogger(__name__)

import demog

def init(transmission_rate=0.065):

    populations = {
        'Population #1' : {
            'n_humans' : 1350,
            'n_infections' : 200,
            'coi' : demog.init_coi_discrete_distribution(n=5, alpha=0.5, beta=8),
            'r0_params': {'r0' : [1, 1],
                          'day' : [0, 25*365]},
            'annual_cycle': lambda t:transmission_rate
            #demog.annual_cycle(0.15,1e-2,8),
                            
        }
    }

    return populations
