import random
import itertools
import numpy as np
import math
import bisect 
import collections

import logging
log = logging.getLogger(__name__)

import utils
import genome as gn
import infection as inf
from human import HumanCohort,HumanIndividual
from migration import MigrationInfo
from scipy.stats import rv_discrete

class Population:
    '''
    An object containing a list of infections and their dynamics
    '''

    def __init__(self, id, parent,
                 n_humans,
                 n_infections=0,
                 migration_rates={},
                 r0_params={'r0' : [1, 1],'day' : [0,20*365]},
                 annual_cycle =  lambda t:0.05
                 ):
                 
        self.id=id
        self.parent=parent
        self.migration_rates = migration_rates
        self.migration_info=MigrationInfo(migration_rates,self.parent.day) # allow for migration start time migration rate at time VERY HACKY....will fix later
        self.susceptibles=HumanCohort(self,n_humans)
        self.infecteds={}
        if n_infections > n_humans:
            raise Exception('Initial infections not to exceed initial humans.')   
        self.annual_cycle=annual_cycle
        self.r0_params = r0_params
        self.vectorial_capacity_fn = self.determine_vectorial_capacity_function(self.r0_params)
        
        for _ in range(n_infections):
            self.add_infection_from_genomes([gn.Genome.from_allele_freq()])
        log.debug(self)

    def __str__(self):
        s  = '%s: '           % self.id
        s += 'humans=%d '     % self.n_humans()
        s += 'infections=%d ' % len(self.infecteds)
        return s        
    
    def add_infection_from_genomes(self,genomes):
        ii=self.add_new_infection(genomes)
        if ii:
            for g in ii.genomes:
                tx=inf.Transmission(infection=ii, genome=g)
                self.notify_transmission([tx])

    def add_new_infection(self,genomes,individual=None):
        if not individual:
            individual=self.susceptibles.pop_individual()
        individual.infection=inf.Infection(individual,genomes)
        self.infecteds[individual.id]=individual
        return individual.infection

    def transmit_infections(self,transmissions):
        n_infections=len(transmissions)
        log.debug('Add %d infections:',n_infections)
        idxs=utils.choose_with_replacement(n_infections,self.n_humans())
        log.debug('Selected individual indices: %s',idxs)
        for idx,transmission in zip(idxs,transmissions):
            genomes=[t.genome for t in transmission]
            if idx<len(self.infecteds):
                i=self.infecteds.values()[idx]
                i.infection.merge_infection(genomes)
                log.debug('Merged strains (idx=%d, id=%d):\n%s',
                          idx,i.id,i.infection)
                self.notify_transmission(transmission,i.infection)
            else:
                log.debug('New infected individual:\n%s',genomes)
                infection=self.add_new_infection(genomes)
                if infection:
                    self.notify_transmission(transmission,infection)

    def notify_transmission(self,transmission,infection=None):
        for t in transmission:
            if infection:
                t.infection=infection
            t.populationId=self.id
            t.day=self.parent.day
        self.parent.notify('infection.transmit',transmission)

    def update(self,dt):
        transmissions=[]
        V=self.vectorial_capacity()
        self.migration_info=MigrationInfo(self.migration_rates,self.parent.day)
        log.info('%s vectorial capacity=%0.2f',self,V)
        for iid,i in self.infecteds.items():
            tt=i.update(dt)
            transmissions.extend(tt)
            if i.infection.expired():
                expired=self.infecteds.pop(iid)
                self.susceptibles.merge_individual(expired)
            elif i.migration.migrating():
                emigrant=self.infecteds.pop(iid)
                self.transmit_emigrant(emigrant)
        if transmissions:
            self.transmit_infections(transmissions)
        self.cohort_migration(dt)
    
    def determine_vectorial_capacity_function(self, r0_params):
        initial_r0 = self.r0_params['r0'][0]
        
        def function(t):
            index =bisect.bisect_left(r0_params['day'], t) 
            base_capacity = initial_r0 * self.annual_cycle(t)
            
            def calculate_vectorial_capacity(start_t, start_R0, end_t, end_R0, t):
                #end_r0 refers to the period in time where R0 stops changing and becomes constant again
                slope = (float(end_R0) - start_R0) / (end_t - start_t)
                delta_t = t-start_t
                # if r0 is constant across the time points, should collapse down to r0 = base_reproduction_rate * r0_scale, where r0_scale = float(current_r0) / initial_r0

                vectorial_capacity = base_capacity * (start_R0 + slope * delta_t) / float(initial_r0)
                return vectorial_capacity
                    
            if index == len(r0_params['day']):
                raise Exception('Cannot extrapolate beyond the given amount of time with which we have a transmission dynamics model for: maximum={max}, tried to do time={t}'.format(max=r0_params['year'][-1], t=t))
            else:
                vectorial_capacity = calculate_vectorial_capacity(r0_params['day'][index - 1], r0_params['r0'][index -1], 
                                    r0_params['day'][index], r0_params['r0'][index],t)
            return vectorial_capacity
        return function
    
    def vectorial_capacity(self):
        return self.vectorial_capacity_fn(t=self.parent.day)

    def n_humans(self):
        return self.susceptibles.n_humans+len(self.infecteds)

    def n_infecteds(self):
        return len(self.infecteds)

    def n_polygenomic(self):
        return sum([i.infection.n_strains()>1 for i in self.infecteds.values()])
    
    def retrieve_infections(self):
        singles, poly = {}, {}
        for i in self.infecteds.values():
            if i.infection.n_strains() == 1:
                singles[i.id] = i.infection.genome_ids()
            elif i.infection.n_strains() > 1:
                poly[i.id] = i.infection.genome_ids()
        return singles, poly
    
    def coi_distribution(self):
        return collections.Counter([i.infection.n_strains() for i in self.infecteds.values()])
        
    def calculate_average_coi(self):
        return np.mean([i.infection.n_strains() for i in self.infecteds.values()])
        
    def transmit_emigrant(self,emigrant):
        src_pop,dest_pop=self.id,emigrant.migration.destination
        self.parent.migrants[dest_pop].append((emigrant,src_pop))

    def cohort_migration(self,dt):
        mig_dest=self.migration_info.destinations_in_timestep
        destinations=mig_dest(self.susceptibles.n_humans,dt)
        self.susceptibles.n_humans -= len(destinations)
        for dest in destinations:
            self.parent.cohort_migrants[dest]+=1
        log.debug('Cohort migration from %s: %s',self.id,destinations)

    def receive_immigrant(self,immigrant,src_pop):
        immigrant.migration=self.migration_info.next_migration()
        # TODO: extend random migration to round-trip concepts using src_pop
        immigrant.parent=self
        self.infecteds[immigrant.id]=immigrant
        
    
