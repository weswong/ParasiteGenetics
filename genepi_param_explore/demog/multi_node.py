import logging
log = logging.getLogger(__name__)

import demog

def init(N=5, V=(0.18,1e-3,8), M=.010):

    def set_population(i):
        if i != 0:
            return { 'n_humans' : 16000,
                     'n_infections' : 1600,
                     'annual_cycle' : lambda t:  0.07, #demog.annual_cycle(*V),
                     'migration_rates' : {'Population #%d'%k:M/N for k in range(N) if k is not i},
                     'coi' : demog.init_coi_discrete_distribution(n=5, alpha=0.5, beta=8),
                     'r0_params': {'r0' : [1, 1, 0.5, 1.1, 1.0],
                                   'day' : [0,5*365,7*365, 8*365, 12*365]
                              } }
        else:
            return { 'n_humans' : 16000,
                     'n_infections' : 1600,
                     'annual_cycle' : lambda t:  0.07, #demog.annual_cycle(*V),
                     'migration_rates' : {'Population #%d'%k:M/N for k in range(N) if k is not i},
                     'coi' : demog.init_coi_discrete_distribution(n=5, alpha=0.5, beta=8),
                     'r0_params': {'r0' : [1, 1, 0.7, 1.0, 1.0],
                                   'day' : [0,5*365,7*365, 8*365, 12*365]
                              } }

    return {'Population #%d'%k:set_population(k) for k in range(N)}
