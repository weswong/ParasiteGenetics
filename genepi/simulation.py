import math
import random
from collections import defaultdict

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import genome as gn
import population as pop

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

    # TODO: Initialization factory - toy topologies, realistic geographies

    populations = {
        'Population #1' : {
            'n_humans' : 500,
            'n_infections' : 20,
            'vectorial_capacity_fn' : annual_cycle(0.2,1e-3,10),
            'migration_rates' : {'Population #2':2e-5},
        },
        'Population #2' : {
            'n_humans' : 500,
            'n_infections' : 150,
            'vectorial_capacity_fn' : annual_cycle(0.045,coeff=0),
            'migration_rates' : {'Population #1':2e-5},
        }
    }

class Params:

    # TODO: Simulation owns instance of parameterization
    #       (e.g. parameter sweeps in parallel)

    working_dir      = 'simulations'
    random_seed      = 8675309

    sim_duration     = 365*10    # days
    sim_tstep        = 21        # days
    incubation       = 25        # days

    @classmethod
    def infectiousness(cls,t):
        if t<cls.incubation:
            return 0
        else:
            mean_prob=0.8*math.exp(-t/50.)+0.05*math.exp(-t/300.)
            return min(1.0,max(0,mean_prob+random.gauss(0,0.1)))

    @classmethod
    def infectious_generator(cls,t=0):
        while True:
            dt=yield cls.infectiousness(t)
            if dt: t+=dt

    @staticmethod
    def get_infection_duration():
        # Maire et al. (2006)
        return random.lognormvariate(5.13,0.8)

    # mutation_prob = ??       # per day, per position

class Simulation:

    def __init__(self, random_seed=Params.random_seed):
        random.seed(random_seed)
        log.debug('Simulation: random seed = %d' % random_seed)
        self.day=0
        self.migrants=defaultdict(list)
        self.cohort_migrants=defaultdict(int)
        self.reports=[]
        self.listeners=defaultdict(list)
        gn.Genome.set_simulation_ref(self)

    def populate_from_demographics(self):
        self.populations={ k:pop.Population(k,self,**v) \
                           for k,v in Demographics.populations.items() }

    def run(self):
        for t in range(Params.sim_duration/Params.sim_tstep):
            self.update()
        for r in self.reports:
            r.write(Params.working_dir)
        for ll in self.listeners.values():
            for l in ll:
                l.write(Params.working_dir)

    def update(self,dt=Params.sim_tstep):
        self.day+=dt
        log.info('\nt=%d'%self.day)
        for p in self.populations.values():
            p.update(dt)
        self.resolve_migration()
        for r in self.reports:
            r.update()

    def resolve_migration(self):
        for dest,emigrant_and_src_list in self.migrants.items():
            for (emigrant,src) in emigrant_and_src_list:
                log.debug('Migrating to %s: infection %s from %s',
                           dest,emigrant.infection,src)
                self.populations[dest].receive_immigrant(emigrant,src)
        self.migrants.clear()
        for dest,n_emigrants in self.cohort_migrants.items():
            self.populations[dest].susceptibles.n_humans+=n_emigrants
        self.cohort_migrants.clear()

    def add_report(self,report_class):
        self.reports.append(report_class(self))

    def add_listener(self,listener_class):
        l=listener_class(self)
        self.listeners[l.event].append(l)

    def notify(self,event,*args):
        listeners=self.listeners.get(event,[])
        for l in listeners:
            l.notify(*args)
