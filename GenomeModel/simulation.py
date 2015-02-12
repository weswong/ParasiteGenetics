import math
import random
import population as pop

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def annual_cycle(year_max, year_min, coeff=2, cycle=365):
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

class Demographics():

    populations = {
        'Test' : {
            'n_humans' : 100,
            'n_infections' : 10,
            'vectorial_capacity' : annual_cycle(0.03,0.005,4),
            #'migration_rates' : {},
        }
    }

class Params():

    sim_duration     = 3*365
    sim_tstep        = 21
    random_seed      = 8675309

    @staticmethod
    def infectiousness(t):
        return 0 if t<25 else 0.05+0.8*math.exp(-t/200.)

    @classmethod
    def infectious_generator(cls,t=0):
        while True:
            dt=yield cls.infectiousness(t)
            if dt: t+=dt

    @staticmethod
    def get_infection_duration():
        return random.uniform(100,300)

    # per-day, per-site synonymous mutation probability
    # mutation_prob = 1e-2

class Simulation():

    day=0

    def __init__(self, random_seed=Params.random_seed):
        random.seed(random_seed)
        log.debug('Simulation: random seed = %d' % random_seed)
        self.populations={ k:pop.Population(k,**v) \
                           for k,v in Demographics.populations.items() }

    def update(self,dt=Params.sim_tstep):
        Simulation.day+=dt
        log.info('t=%d'%self.day)
        for p in self.populations.values():
            p.update(dt)

    def run(self):
        for t in range(Params.sim_duration/Params.sim_tstep):
            self.update()
