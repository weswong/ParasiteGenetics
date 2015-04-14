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
            'vectorial_capacity_fn' : sim.annual_cycle(0.2,2e-3,10)
        }
    }

    sim.Demographics.populations=populations
