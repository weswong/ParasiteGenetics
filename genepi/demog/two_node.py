import os

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import genepi.simulation as sim

def init():

    populations = {
        'Population #1' : {
            'n_humans' : 500,
            'n_infections' : 20,
            'vectorial_capacity_fn' : sim.annual_cycle(0.2,1e-3,10),
            'migration_rates' : {'Population #2':2e-5},
        },
        'Population #2' : {
            'n_humans' : 500,
            'n_infections' : 150,
            'vectorial_capacity_fn' : sim.annual_cycle(0.045,coeff=0),
            'migration_rates' : {'Population #1':2e-5},
        }
    }

    sim.Demographics.populations=populations
