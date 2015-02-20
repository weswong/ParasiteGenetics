import math
import random
import itertools
import genome as gn
import simulation as sim

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

max_transmit_strains=10

def sample_n_oocysts():
    # Fit A.Ouedraogo membrane feeding data from Burkina Faso
    # Private communication, in preparation (2015)
    return 1 + int(random.weibullvariate(alpha=2.5,beta=0.65))

def sample_n_hepatocytes():
    # Bejon et al. "Calculation of Liver-to-Blood Inocula..." (2005)
    # ln(5)~=1.6, ln(2.7)~=1
    return max(1,int(random.lognormvariate(mu=1.6,sigma=0.8)))

def accumulate_cdf(iterable):
    cdf,subtotal=[],0
    norm=float(sum(iterable))
    if not norm:
        return []
    for it in iterable:
        subtotal += it/norm
        cdf.append(subtotal)
    cdf[-1]=1.0
    return cdf

def nextTime(rateParameter):
    if rateParameter<=0:
        raise Exception('Cannot calculate next time from zero rate.')
    return -math.log(1.0 - random.random()) / rateParameter

def poissonRandom(rateParameter):
    if rateParameter<=0:
        return 0
    sumTime=0
    N=0
    while True:
        sumTime+=nextTime(rateParameter)
        if sumTime>1:
            break
        N+=1
    return N

def weighted_choice(cumwts):
    R = random.random()
    idx=sum(itertools.takewhile(bool, (cw < R for cw in cumwts)))
    return idx

class MigratingIndividual:
    def __init__(self,in_days=[],destination=[]):
        self.in_days=in_days
        self.destination=destination

    def __str__(self):
        s  = 'destination=%s: ' % self.destination
        s += 'in_days=%f '      % self.in_days
        return s

    def update(self,dt):
        if self.in_days:
            self.in_days -= dt

    def migrating(self):
        return True if self.in_days and self.in_days<=0 else False

class Infection:
    '''
    An infected individual with one or more parasite strains,
    including the dynamics of parasite propagation through
    the mosquito vector into a new human host.
    '''

    id=itertools.count()

    def __init__(self,genomes=[]):
        self.id=Infection.id.next()
        log.debug('Infection: id=%d', self.id)
        self.set_infection_timers()
        self.genomes=genomes
        log.debug('%s', self)
        self.migration=MigratingIndividual()

    def __repr__(self):
        return str(self.id)

    def __str__(self):
        return '\n'.join([str(g) for g in self.genomes])

    @classmethod
    def from_random(cls,n_clones=1):
        return cls([gn.Genome.from_allele_frequencies() for _ in range(n_clones)])

    def set_infection_timers(self,t=0):
        self.infectiousness=sim.Params.infectious_generator(t)
        self.infectiousness.send(None)
        self.infection_timer=sim.Params.get_infection_duration()

    def update(self,dt,vectorial_capacity):
        self.infection_timer  -= dt
        log.debug('infection_timer=%d',self.infection_timer)
        self.migration.update(dt)
        self.infectiousness.send(dt)
        transmit_rate = vectorial_capacity*dt*next(self.infectiousness)
        n_transmit=poissonRandom(transmit_rate)
        log.debug('transmit_rate=%0.2f  n_transmit=%d',transmit_rate,n_transmit)
        transmits=[self.transmit() for _ in range(n_transmit)]
        return transmits

    def transmit(self):
        n_hep,n_ooc=sample_n_hepatocytes(),sample_n_oocysts()
        log.debug('Sample %d hepatocyte(s) from %d oocyst(s):',n_hep,n_ooc)
        if n_hep > max_transmit_strains:
            log.debug('Truncating to %d hepatocytes:',max_transmit_strains)
            n_hep=max_transmit_strains
        sporozoite_strains=[]
        for g1,g2 in self.sample_gametocyte_pairs(n_ooc):
            meiotic_products=gn.meiosis(g1,g2)
            log.debug([str(mp) for mp in meiotic_products])
            sporozoite_strains.extend(meiotic_products)
        hepatocytes=[random.choice(sporozoite_strains) for _ in range(n_hep)]
        return Infection(gn.distinct(hepatocytes))

    def sample_gametocyte_pairs(self, N):
        pairs=[]
        w=self.gametocyte_strain_cdf()
        for _ in range(N):
            pairs.append([self.select_strain(w),self.select_strain(w)])
        return pairs

    def select_strain(self,cumwts):
        return self.genomes[weighted_choice(cumwts)]

    def gametocyte_strain_cdf(self):
        # TODO: something more skewed
        #       to account for blood-stage dynamics
        return accumulate_cdf([1]*self.n_strains())

    def n_strains(self):
        return len(self.genomes)

    def add_infection(self,inf):
        self.genomes.extend(inf.genomes)
        self.genomes=gn.distinct(self.genomes)
        self.set_infection_timers(t=sim.Params.incubation)
