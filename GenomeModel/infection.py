import itertools
import random
import genome as gn
from utils import log

def sample_n_oocysts():
    # Fit A.Ouedraogo membrane feeding data from Burkina Faso
    # Private communication, in preparation (2015)
    return 1 + int(random.weibullvariate(alpha=2.5,beta=0.65))

def sample_n_hepatocytes():
    # Bejon et al (2005)
    # ln(5)~=1.6, ln(2.7)~=1
    return int(random.lognormvariate(mu=1.6,sigma=1))

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

    def transmit(self):
        pass

    def sample_gametocyte_pairs(self, N):
        pairs=[]
        w=self.gametocyte_strain_cdf()
        for _ in range(N):
            pairs.append([self.weighted_choice(w),self.weighted_choice(w)])
        return pairs

    def weighted_choice(self, cumwts):
        R = random.random()
        idx=sum(itertools.takewhile(bool, (cw < R for cw in cumwts)))
        return self.genomes[idx]

    def gametocyte_strain_cdf(self):
        # TODO: something more skewed 
        #       to account for blood-stage dynamics
        total=0
        cdf=[]
        n_strains=self.n_strains()
        for _ in range(n_strains):
            total += 1.0/n_strains
            cdf.append(total)
        cdf[-1]=1.0
        return cdf

    def n_strains(self):
        return len(self.genomes)
