import math
import random
import itertools
from collections import defaultdict

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import utils
import genome as gn
from human import HumanIndividual
import simulation as sim

max_transmit_strains=10

def sample_n_oocysts():
    # Fit A.Ouedraogo membrane feeding data from Burkina Faso
    # Private communication, in preparation (2015)
    return 1 + int(random.weibullvariate(alpha=2.5,beta=0.65))

def sample_n_hepatocytes():
    # Bejon et al. "Calculation of Liver-to-Blood Inocula..." (2005)
    # ln(5)~=1.6, ln(2.7)~=1
    return max(1,int(random.lognormvariate(mu=1.6,sigma=0.8)))

def sample_oocyst_products(n_hep,n_ooc):
    oocyst_product_idxs=list(itertools.product(range(n_ooc),range(4)))
    hep_idxs=[random.choice(oocyst_product_idxs) for _ in range(n_hep)]
    product_idxs=defaultdict(set)
    for o_idx,m_idx in hep_idxs:
        product_idxs[o_idx].add(m_idx)
    n_products_by_oocyst = [len(v) for k,v in product_idxs.items()]
    log.debug('meiotic products to be sampled per oocyst:%s'%n_products_by_oocyst)
    return n_products_by_oocyst

class Transmission:
    '''
    Characteristics of transmitted sporozoites and parent gametocytes
    '''
    def __init__(self,parentGenomeIds=(None,None),genome=None,
                      parentInfectionId=None,infectionId=None,
                      populationId=None,day=None):
        # sporozoite properties
        self.genome=genome
        self.infectionId=infectionId
        # male + female gametocyte properties
        self.parentGenomeIds=parentGenomeIds
        self.parentInfectionId=parentInfectionId
        # other info
        self.populationId=populationId
        self.day=day

    def to_tuple(self):
        return (self.day,self.populationId,
                self.infectionId,self.parentInfectionId,
                self.parentGenomeIds[0],self.parentGenomeIds[1],
                self.genome.id)

class Infection():
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
        self.infectiousness=sim.Params.infectious_generator(t)
        self.infectiousness.send(None)
        self.infection_timer=sim.Params.get_infection_duration()

    def update(self,dt,vectorial_capacity):
        self.infection_timer  -= dt
        self.infectiousness.send(dt)
        transmit_rate = vectorial_capacity*dt*next(self.infectiousness)
        n_transmit=utils.poissonRandom(transmit_rate)
        log.debug('  id=%d: infection_timer=%d  transmit_rate=%0.2f  n_transmit=%d',
                  self.id,self.infection_timer,transmit_rate,n_transmit)
        transmissions=[self.transmit() for _ in range(n_transmit)]
        return transmissions

    def transmit(self):
        if self.n_strains() == 1:
            clone=self.genomes[0]
            log.debug('Clonal transmission of genome id=%d'%clone.id)
            return [Transmission((clone.id,clone.id),clone,self.id)]
        n_hep,n_ooc=sample_n_hepatocytes(),sample_n_oocysts()
        log.debug('Sample %d hepatocyte(s) from %d oocyst(s):',n_hep,n_ooc)
        if n_hep > max_transmit_strains:
            log.debug('Truncating to %d hepatocytes:',max_transmit_strains)
            n_hep=max_transmit_strains
        n_products=sample_oocyst_products(n_hep,n_ooc)
        gametocyte_pairs=self.sample_gametocyte_pairs(len(n_products))
        sporozoites=gn.distinct_sporozoites_from(gametocyte_pairs,n_products)
        for s in sporozoites:
            s.parentInfectionId=self.id
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
        #       to account for blood-stage dynamics
        return utils.accumulate_cdf([1]*self.n_strains())

    def n_strains(self):
        return len(self.genomes)

    def individual(self):
        return self.parent

    def population(self):
        return self.parent.parent

    def simulation(self):
        return self.parent.parent.parent

    def add_infection(self,genomes):
        self.genomes.extend(genomes)
        self.genomes=gn.distinct(self.genomes)
        self.set_infection_timers(t=sim.Params.incubation)
