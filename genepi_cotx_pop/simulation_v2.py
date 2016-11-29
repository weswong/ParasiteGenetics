import math
import random
from Queue import PriorityQueue
from collections import defaultdict
from importlib import import_module

import logging
log = logging.getLogger(__name__)
print log

import genome as gn
import infection as inf
import population as pop
import sys

class Simulation:

    class Params:
        def __init__(self,**kwargs):
            self.working_dir  = 'simulations'
            self.random_seed  = 8675309

            self.sim_duration = 365*5    # days
            self.sim_tstep    = 21        # days

            #infection 
            self.n_oocyst_params = (2.5, 0.65)     #weibull distribution, (lambda, shape_parameter)
            self.n_hepatocyte_params = (1.6, 0.8)  #lognormal distribution, (mu*, stdev*)
            self.infection_duration_params = (5.13, 0.8) #lognormal distribution, (mu, stdev)
            
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
        self.events=PriorityQueue()
        inf.Infection.sample_n_oocysts_function(*self.params.n_oocyst_params)
        inf.Infection.sample_n_hepatocytes_function(*self.params.n_hepatocyte_params)
        inf.Infection.get_infection_duration_function(*self.params.infection_duration_params)
        gn.Genome.set_simulation_ref(self)
    
    
    def populate_from_demographics(self,mod='single_node',*args,**kwargs):
        try:
            mod=import_module('.'.join(['genepi','demog',mod]))
            demog=mod.init(*args,**kwargs)
        except ImportError as e:
            sys.exit("ImportError for demog_source: %s"%e)
        self.populations={ k:pop.Population(k,self,**v) \
                           for k,v in demog.items() }
        self.demog_params=demog.items()

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
        while True:
            if self.events.empty() or self.events.queue[0][0] > self.day:
                break
            evt = self.events.get()[1]
            log.info('Executing event.')
            evt(self)
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

    def add_event(self,day,event):
        self.events.put((day,event))

    def add_reports(self,*args):
        for report_class in args:
            self.reports.append(report_class(self))

    def add_listeners(self,*args):
        for listener_class in args:
            l=listener_class(self)
            self.listeners[l.event].append(l)

    def notify(self,event,*args):
        listeners=self.listeners.get(event,[])
        for l in listeners:
            l.notify(*args)
