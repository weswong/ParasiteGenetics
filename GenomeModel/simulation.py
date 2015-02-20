import math
import random
from collections import defaultdict
import population as pop

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def annual_cycle(year_max, year_min=0, coeff=2, cycle=365):
    def f(t):
        if coeff == 0:
            return year_max
        if coeff % 2:
            raise Exception('Only supporting even sinusoid powers.')
        return year_min+(year_max-year_min)*pow(math.cos(3.1416*t/cycle),coeff)
    return f

def gravity(p1,p2,d12,G=1e-3):
    # (1e-3) : 1/day to pop=1,000 village @ 1km
    return G*p1*p2/d12**2

class Demographics:

    # TODO: initialization factory

    populations = {
        'Test1' : {
            'n_humans' : 100,
            'n_infections' : 0,
            'vectorial_capacity' : annual_cycle(0.03,0.005,4),
            'migration_rates' : {'Test2':0.1},
        },
        'Test2' : {
            'n_humans' : 100,
            'n_infections' : 10,
            'vectorial_capacity' : annual_cycle(0.01,coeff=0),
            'migration_rates' : {'Test1':0.1},
        }
    }

class Params:

    sim_duration     = 365       # days
    sim_tstep        = 21        # days
    incubation       = 25        # days
    random_seed      = 8675309

    @classmethod
    def infectiousness(cls,t):
        return 0 if t<cls.incubation else 0.05+0.8*math.exp(-t/200.)

    @classmethod
    def infectious_generator(cls,t=0):
        while True:
            dt=yield cls.infectiousness(t)
            if dt: t+=dt

    @staticmethod
    def get_infection_duration():
        return random.uniform(100,300)

    # mutation_prob = ??       # per day, per position

class Simulation:

    def __init__(self, random_seed=Params.random_seed):
        random.seed(random_seed)
        log.debug('Simulation: random seed = %d' % random_seed)
        self.day=0
        self.populations={ k:pop.Population(k,self,**v) \
                           for k,v in Demographics.populations.items() }
        self.migrants=defaultdict(list)

    def run(self):
        for t in range(Params.sim_duration/Params.sim_tstep):
            self.update()

    def update(self,dt=Params.sim_tstep):
        self.day+=dt
        log.info('t=%d'%self.day)
        for p in self.populations.values():
            p.update(dt)
        self.resolve_migration()

    def resolve_migration(self):
        for dest,emigrant_and_src_list in self.migrants.items():
            for emigrant_and_src in emigrant_and_src_list:
                log.debug('Migrating to %s: infection %s from %s',dest,*emigrant_and_src)
                self.populations[dest].receive_immigrant(*emigrant_and_src)
        self.migrants.clear()
