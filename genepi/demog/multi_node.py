import os

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import genepi.simulation as sim

def init(N=10, V=(0.18,1e-3,8), M=5e-4):

    def set_population(i):
        return { 'n_humans' : 500,
                 'n_infections' : 20,
                 'vectorial_capacity_fn' : sim.annual_cycle(*V),
                 'migration_rates' : {'Population #%d'%k:M/N for k in range(N) if k is not i} }

    populations = {'Population #%d'%k:set_population(k) for k in range(N)}

    sim.Demographics.populations=populations
