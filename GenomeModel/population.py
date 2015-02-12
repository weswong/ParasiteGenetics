import random
import itertools
import genome as gn
import infection as inf
import simulation as sim

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def choose_without_replacement(M,N):
    '''
    O(M) in choose M from N scenario,
    which is much faster for typical use case
    than random.sample, which is O(N)
    '''
    if M==N:
        return range(M)
    chosen_idxs=set()
    for j in range(N-M,N):
        t = random.randint(0,j-1)
        idx = t if t not in chosen_idxs else j
        chosen_idxs.add(idx)
    return list(chosen_idxs)

class Population():
    '''
    An object containing a list of infections and their dynamics
    '''

    def __init__(self, id,
                 n_humans,
                 n_infections=0,
                 vectorial_capacity=lambda t:0.05,
                 migration_rates={}):
        self.id=id
        self.n_humans=n_humans
        self.infections=[inf.Infection.from_random() for _ in range(n_infections)]
        self.vectorial_capacity=vectorial_capacity
        self.migration_rates=migration_rates
        log.debug(self)

    def __repr__(self):
        s  = '%s: ' % self.id
        s += 'humans=%d ' % self.n_humans
        s += 'infections=%d ' % len(self.infections)
        return s

    def add_infections(self,infections):
        n_infections=len(infections)
        log.debug('Add %d infections:',n_infections)
        log.debug('\n\n'.join([str(i) for i in infections]))
        idxs=choose_without_replacement(n_infections,self.n_humans)
        log.debug('Selected individual indices: %s',idxs)
        for idx,infection in zip(idxs,infections):
            if idx<len(self.infections):
                self.infections[idx].add_infection(infection)
                log.debug('Merged strains (idx=%d):\n%s',idx,self.infections[idx])
            else:
                log.debug('New infected individual:\n%s',infection)
                self.infections.append(infection)

    def update(self,dt):
        transmissions=[]
        V=self.vectorial_capacity(sim.Simulation.day)
        #log.debug('Vectorial capacity = %0.2f',V)
        for i in self.infections:
            transmit=i.update(dt,V)
            if transmit:
                transmissions.append(transmit)
        if transmissions:
            self.add_infections(transmissions)
        self.infections[:] = [i for i in self.infections if i.infection_timer>0]
        log.debug('%s vectorial capacity=%0.2f',self,V)
