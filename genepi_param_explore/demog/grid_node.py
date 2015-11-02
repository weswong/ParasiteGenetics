from itertools import product

import logging
log = logging.getLogger(__name__)

import demog

def init(L=3, V=(0.18,1e-3,8), M=5e-4):

    def idx(i,j):
        return i*L+j

    def in_range(x):
        return 0 <= x < L

    def not_same(di, dj):
        return not (di==0 and dj==0)

    def set_population(i,j):
        dx = [-1, 0, 1]
        return { 'n_humans' : 500,
                 'n_infections' : 20,
                 'annual_cycle' : lambda t: 0.0257 ＃demog.annual_cycle(*V),
                 'migration_rates' : {'Population #%d'%idx(i+di,j+dj):M/8 for (di,dj) in product(dx,dx)
                                      if not_same(di,dj) and in_range(i+di) and in_range(j+dj)} }

    return {'Population #%d'%idx(i,j):set_population(i,j) for (i,j) in product(range(L),range(L))}
