import itertools
import genome as gn
from utils import log

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
        self.genomes=genomes
        log.debug('%s', self)

    def __repr__(self):
        return '\n'.join([str(g) for g in self.genomes])
