import random
import itertools

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import utils
import genome as gn
import infection as inf
from migration import MigrationInfo

class Population:
    '''
    An object containing a list of infections and their dynamics
    '''

    def __init__(self, id, parent,
                 n_humans,
                 n_infections=0,
                 vectorial_capacity=lambda t:0.05,
                 migration_rates={}):
        self.id=id
        self.parent=parent
        self.vectorial_capacity=vectorial_capacity
        self.migration_info=MigrationInfo(migration_rates)
        self.n_humans=n_humans
        self.infections={}
        for _ in range(n_infections):
            i=inf.Infection.from_random()
            self.add_new_infection(i)
        log.debug(self)

    def __str__(self):
        s  = '%s: '           % self.id
        s += 'humans=%d '     % self.n_humans
        s += 'infections=%d ' % len(self.infections)
        return s

    def add_new_infection(self,infection):
        infection.migration=self.migration_info.next_migration()
        self.infections[infection.id]=infection

    def add_infections(self,infections):
        n_infections=len(infections)
        log.debug('Add %d infections:',n_infections)
        log.debug('\n\n'.join([str(i) for i in infections]))
        idxs=utils.choose_with_replacement(n_infections,self.n_humans)
        log.debug('Selected individual indices: %s',idxs)
        for idx,infection in zip(idxs,infections):
            if idx<len(self.infections):
                i=self.infections.values()[idx]
                i.add_infection(infection)
                log.debug('Merged strains (idx=%d, id=%d):\n%s',idx,i.id,i)
            else:
                log.debug('New infected individual:\n%s',infection)
                self.add_new_infection(infection)

    def update(self,dt):
        transmissions=[]
        V=self.vectorial_capacity(self.parent.day)
        for iid,i in self.infections.items():
            transmits=i.update(dt,V)
            transmissions.extend(transmits)
            if i.infection_timer<=0:
                expired=self.infections.pop(iid)
            elif i.migration.in_days<=0:
                emigrant=self.infections.pop(iid)
                self.transmit_emigrant(emigrant)
        if transmissions:
            self.add_infections(transmissions)
        log.debug('%s vectorial capacity=%0.2f',self,V)
        # N.B. this line doesn't report currently emigrating
        #      who are held in limbo by the Simulation object
        #      until being redistributed at end of update()

    def transmit_emigrant(self,emigrant):
        src_pop,dest_pop=self.id,emigrant.migration.destination
        self.parent.migrants[dest_pop].append((emigrant,src_pop))
        self.n_humans-=1

    def receive_immigrant(self,immigrant,src_pop):
        immigrant.migration=self.migration_info.next_migration()
        # TODO: extend random migration to round-trip concepts using src_pop
        self.infections[immigrant.id]=immigrant
        self.n_humans+=1
