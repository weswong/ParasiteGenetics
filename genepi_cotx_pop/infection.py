import math
import random
import itertools
from collections import defaultdict

import logging
log = logging.getLogger(__name__)

import numpy as np

import utils
import genome as gn
from human import HumanIndividual
from transmission import Transmission

max_transmit_strains=10

incubation=25 # days

def infectiousness(t):
    if t<incubation:
        return 0
    else:
        # TODO: choose functional form based on EMOD DTK calibration
        mean_prob = 0.8 * math.exp(-t/50.) + 0.05 * math.exp(-t/300.)
        return min(1.0, max(1e-6, mean_prob+random.gauss(0, 0.1)))

def infectious_generator(t=0):
    while True:
        dt=yield infectiousness(t)
        if dt: t+=dt

def sample_oocyst_products(n_hep,n_ooc):
    oocyst_product_idxs=list(itertools.product(range(n_ooc),range(4)))
    hep_idxs=[random.choice(oocyst_product_idxs) for _ in range(n_hep)]
    product_idxs=defaultdict(set)
    for o_idx,m_idx in hep_idxs:
        product_idxs[o_idx].add(m_idx)
    n_products_by_oocyst = [len(v) for k,v in product_idxs.items()]
    log.debug('meiotic products to be sampled per oocyst:%s'%n_products_by_oocyst)
    return n_products_by_oocyst
    
class Infection:
    '''
    An infection in a human containing one or more parasite strains,
    including the dynamics of parasite propagation through
    the mosquito vector into a new human host.
    '''

    id=itertools.count()

    def __init__(self,parent,genomes=[]):
        self.id=Infection.id.next()
        log.debug('Infection: id=%d', self.id)
        self.parent=parent
        self.set_infection_timers()
        self.genomes=genomes
        log.debug('%s', self)

    def __str__(self):
        return '\n'.join([str(g) for g in self.genomes])

    def set_infection_timers(self,t=0):
        self.infectiousness=infectious_generator(t)
        self.infectiousness.send(None)
        self.infection_timer=self.infection_duration()

    def update(self,dt,vectorial_capacity):
        self.infection_timer  -= dt
        self.infectiousness.send(dt)
        fitness = [g.fitness() for g in self.genomes]
        if not fitness:
            return [] # no genomes (e.g. from drug clearance)
        mean_fitness = np.mean(fitness) # TODO: something other than mean?
        transmit_rate = vectorial_capacity*dt*next(self.infectiousness)*mean_fitness
        n_transmit=utils.poissonRandom(transmit_rate)
        log.debug('  id=%d: infection_timer=%d  transmit_rate=%0.2f  n_transmit=%d',
                  self.id,self.infection_timer,transmit_rate,n_transmit)
        transmissions=[self.transmit() for _ in range(n_transmit)]
        return transmissions

    def transmit(self):
        if self.n_strains() == 1:
            clone=self.genomes[0]
            log.debug('Clonal transmission of genome id=%d'%clone.id)
            return [Transmission((clone.id,clone.id),clone,self)]
        n_hep,n_ooc=self.sample_n_hep(), self.sample_n_oocysts()
        log.debug('Sample %d hepatocyte(s) from %d oocyst(s):',n_hep,n_ooc)
        if n_hep > max_transmit_strains:
            log.debug('Truncating to %d hepatocytes:',max_transmit_strains)
            n_hep=max_transmit_strains
        n_products=sample_oocyst_products(n_hep,n_ooc)
        gametocyte_pairs=self.sample_gametocyte_pairs(len(n_products))
        sporozoites=gn.distinct_sporozoites_from(gametocyte_pairs,n_products)
        for s in sporozoites:
            s.parentInfection=self
        return sporozoites

    def sample_gametocyte_pairs(self, N):
        pairs=[]
        w=self.gametocyte_strain_cdf()
        for _ in range(N):
            pairs.append([self.select_strain(w),self.select_strain(w)])
        return pairs

    def select_strain(self,cumwts):
        return self.genomes[utils.weighted_choice(cumwts)]

    def expired(self):
        return self.infection_timer<=0

    def gametocyte_strain_cdf(self):
        # TODO: something more skewed
        #       to account for blood-stage dynamics, e.g. in EMOD DTK
        #return utils.accumulate_cdf([1]*self.n_strains()) # random
        return utils.accumulate_cdf([g.fitness() for g in self.genomes]) # fitness weighted

    def n_strains(self):
        return len(self.genomes)

    def merge_infection(self,genomes):
        self.genomes.extend(genomes)
        self.genomes=gn.distinct(self.genomes)
        self.set_infection_timers(t=incubation)
    
    def genome_ids(self):
        return [genome.id for genome in self.genomes]
    
    @classmethod
    def sample_n_oocysts_function(cls, alpha=2.5,beta=1.0):
         # Fit A.Ouedraogo membrane feeding data from Burkina Faso
        # Private communication, in preparation (2015)
        def f(self):
            return 1 + int(random.weibullvariate(alpha, beta))
        cls.sample_n_oocysts = f
    
    @classmethod
    def sample_n_hepatocytes_function(cls, mu=1.8,sigma=0.8):
        # Bejon et al. "Calculation of Liver-to-Blood Inocula..." (2005)
        # ln(5)~=1.6, ln(2.7)~=1
        def f(self):
            return max(1,int(random.lognormvariate(mu,sigma)))
        cls.sample_n_hep = f
    
    @classmethod
    def get_infection_duration_function(cls, mu=5.13, sigma=0.8):  # Maire et al. (2006)
        def f(self):
            return random.lognormvariate(mu, sigma)
        cls.infection_duration = f
    #return random.lognormvariate(3.2, 1.55)
