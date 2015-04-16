import os
import logging
log = logging.getLogger(__name__)

import demog

def init(N=10, V=(0.18,1e-3,8), M=5e-4):

    def set_population(i):
        return { 'n_humans' : 500,
                 'n_infections' : 20,
                 'vectorial_capacity_fn' : demog.annual_cycle(*V),
                 'migration_rates' : {'Population #%d'%k:M/N for k in range(N) if k is not i} }

    return {'Population #%d'%k:set_population(k) for k in range(N)}
