import math
import random
from collections import defaultdict
from importlib import import_module

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import genome as gn
import population as pop

class Simulation:

    class Params:
        def __init__(self,**kwargs):
            self.working_dir  = 'simulations'
            self.random_seed  = 8675309

            self.sim_duration = 365*10    # days
            self.sim_tstep    = 21        # days

            for k,v in kwargs.items():
                try:
                    d=getattr(self,k)
                    setattr(self,k,v)
                    log.info('  Setting %s to %s'%(k,v))
                except AttributeError:
                    print('Cannot override non-existent parameter: %s' % k)

    def __init__(self, **kwargs):
        self.params=self.Params(**kwargs)
        random.seed(self.params.random_seed)
        log.debug('Simulation: random seed = %d' % self.params.random_seed)
        self.day=0
        self.migrants=defaultdict(list)
        self.cohort_migrants=defaultdict(int)
        self.reports=[]
        self.listeners=defaultdict(list)
        gn.Genome.set_simulation_ref(self)

    def populate_from_demographics(self,mod='single_node',*args,**kwargs):
        try:
            mod=import_module('.'.join(['genepi','demog',mod]))
            demog=mod.init(*args,**kwargs)
        except ImportError as e:
            sys.exit("ImportError for demog_source: %s"%e)
        self.populations={ k:pop.Population(k,self,**v) \
                           for k,v in demog.items() }

    def run(self):
        for t in range(self.params.sim_duration/self.params.sim_tstep):
            self.update()
        for r in self.reports:
            r.write(self.params.working_dir)
        for ll in self.listeners.values():
            for l in ll:
                l.write(self.params.working_dir)

    def update(self):
        dt=self.params.sim_tstep
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
