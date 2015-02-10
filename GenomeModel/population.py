import random
import itertools
import genome as gn
import infection as inf

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

class Population():
    '''
    An object containing a list of infections and their dynamics
    '''

    def __init__(self, id, 
                 n_humans=0, 
                 n_infections=0, 
                 migration_rates={}):
        self.id=id
        self.n_humans=n_humans
        self.migration_rates=migration_rates
        self.infections=[inf.Infection.from_random() for _ in range(n_infections)]
        log.debug(self)

    def __repr__(self):
        s='Node %d:' % self.id
        s+='\thumans=%d' % self.n_humans
        s+='\tinfections=%d' % len(self.infections)
        return s

    def add_infections(self,n_infections):
        log.debug('Add %d infections',n_infections)
        for _ in range(n_infections):
            # decide whether new infection or reinfection
            # if new:
            #   add
            # else:
            #   pick existing infection
            #   merge
            pass

    def update(self):
        pass