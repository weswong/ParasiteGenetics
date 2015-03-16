import random
import itertools

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import utils
import genome as gn
import infection as inf
from human import HumanCohort,HumanIndividual
from migration import MigrationInfo

class Population:
    '''
    An object containing a list of infections and their dynamics
    '''

    def __init__(self, id, parent,
                 n_humans,
                 n_infections=0,
                 vectorial_capacity_fn=lambda t:0.05,
                 migration_rates={}):
        self.id=id
        self.parent=parent
        self.vectorial_capacity_fn=vectorial_capacity_fn
        self.migration_info=MigrationInfo(migration_rates)
        self.susceptibles=HumanCohort(self,n_humans)
        self.infecteds={}
        if n_infections > n_humans:
            raise Exception('Initial infections not to exceed initial humans.')
        for _ in range(n_infections):
            i=inf.Infection.from_random()
            self.add_new_infection(i)
        log.debug(self)

    def __str__(self):
        s  = '%s: '           % self.id
        s += 'humans=%d '     % self.n_humans()
        s += 'infections=%d ' % len(self.infecteds)
        return s

    def add_new_infection(self,infection,individual=None):
        if not individual:
            individual=self.susceptibles.pop_individual()
        individual.infection=infection
        self.infecteds[individual.id]=individual

    def add_infections(self,infections):
        n_infections=len(infections)
        log.debug('Add %d infections:',n_infections)
        #log.debug('\n\n'.join([str(i) for i in infections]))
        idxs=utils.choose_with_replacement(n_infections,self.n_humans())
        log.debug('Selected individual indices: %s',idxs)
        for idx,infection in zip(idxs,infections):
            if idx<len(self.infecteds):
                i=self.infecteds.values()[idx]
                i.infection.add_infection(infection)
                log.debug('Merged strains (idx=%d, id=%d):\n%s',
                          idx,i.id,i.infection)
            else:
                log.debug('New infected individual:\n%s',infection)
                self.add_new_infection(infection)

    #@profile
    def update(self,dt):
        transmissions=[]
        V=self.vectorial_capacity()
        log.info('%s vectorial capacity=%0.2f',self,V)
        for iid,i in self.infecteds.items():
            transmits=i.update(dt)
            transmissions.extend(transmits)
            if i.infection.expired():
                expired=self.infecteds.pop(iid)
                self.susceptibles.merge_individual(expired)
            elif i.migration.migrating():
                emigrant=self.infecteds.pop(iid)
                self.transmit_emigrant(emigrant)
        if transmissions:
            self.add_infections(transmissions)
        self.cohort_migration(dt)

    def vectorial_capacity(self):
        return self.vectorial_capacity_fn(self.parent.day)

    def n_humans(self):
        return self.susceptibles.n_humans+len(self.infecteds)

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
