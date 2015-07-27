import logging
log = logging.getLogger(__name__)

import demog

def init():

    populations = {
        'Population #1' : {
            'n_humans' : 500,
            'n_infections' : 20,
            'vectorial_capacity_fn' : demog.annual_cycle(0.2,2e-3,10)
            'coi' : init_coi_discrete_distribution()
        }
    }

    return populations
