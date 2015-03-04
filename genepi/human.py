import itertools

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

from migration import Migration

class HumanCohort:
    '''
    A cohort used to track aggregate properties of uninfected humans
    '''
    def __init__(self,parent,n_humans):
        self.parent=parent
        self.n_humans=n_humans

    def merge_individual(self,individual):
        self.n_humans+=1

    def pop_individual(self):
        self.n_humans-=1
        return HumanIndividual(self.parent)

class HumanIndividual:
    '''
    An individual instance for tracking infected humans
    '''
    id=itertools.count()

    def __init__(self,parent,infection=None):
        self.id=HumanIndividual.id.next()
        log.debug('HumanIndividual: id=%d', self.id)
        self.parent=parent
        self.infection=infection
        self.migration=parent.migration_info.next_migration()

    def update(self,dt):
        self.migration.update(dt)
        return self.infection.update(dt,self.parent.vectorial_capacity())
