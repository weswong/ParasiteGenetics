import random
import itertools
import genome as gn
import infection as inf

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def choose_with_replacement(M,N):
    choose = random.choice
    indices = range(N)
    return [choose(indices) for _ in range(M)]

def choose_without_replacement(M,N):
    '''
    O(M) in choose M from N scenario,
    which is much faster for typical use case
    than random.sample, which is O(N)
    '''
    if M>N:
        raise Exception('Cannot sample %d from %d without replacement',(M,N))
    if M==N:
        return range(M)
    chosen_idxs=set()
    for j in range(N-M,N):
        t = random.randint(0,j)
        idx = t if t not in chosen_idxs else j
        chosen_idxs.add(idx)
    return list(chosen_idxs)

class MigrationInfo:
    def __init__(self,rates):
        self.destinations=rates.keys()
        self.relative_rates=inf.accumulate_cdf(rates.values())
        self.total_rate=sum(rates.values())

    def __str__(self):
        s  = 'destinations=%s: '  % self.destinations
        s += 'relative rates=%s ' % self.relative_rates
        s += 'total_rate=%f '     % self.total_rate
        return s

    def next_migration(self):
        if not self.total_rate:
            return inf.MigratingIndividual() # nowhere to migrate
        in_days=random.expovariate(self.total_rate)
        destination=self.destinations[inf.weighted_choice(self.relative_rates)]
        return inf.MigratingIndividual(in_days, destination)

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

    def __repr__(self):
        return str(self.id)

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
        idxs=choose_with_replacement(n_infections,self.n_humans)
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
