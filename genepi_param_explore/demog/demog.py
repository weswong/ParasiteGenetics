import math
import itertools
import bisect
from scipy.special import gamma 
from scipy.stats import beta as beta_distribution

def annual_cycle(year_max, year_min=0, coeff=2, cycle=365):
    def f(t):
        if coeff == 0:
            return year_max
        if coeff % 2:
            raise Exception('Only supporting even sinusoid powers.')
        return year_min+(year_max-year_min)*pow(math.cos(3.1416*t/cycle),coeff)
    return f
    

def gravity(p1,p2,d12,G=1e-3):
    # (1e-3) : 1/day to pop=1,000 village @ 1km
    return G*p1*p2/d12**2

def init_coi_discrete_distribution(n=5, alpha=0.5, beta=8):
    '''feed in an expectation for the distribution of COI in a given population. 
    using a bastardized version of the beta_binomial density... PLACEHOLDER
    distribution eyeballed off of Galinksy et al 2015, COIL paper
    success = mix, failure = not-mixed, so k=0 refers to single-infection, k=1 refers to COI = 2, etc etc (necessarily assumes that COI estimates are independent of one another)
   '''
    
    coi_discrete_distribution = {}
    for k in range(n):
        # dictionary is automatially keyed to add 1 to k to make it more intuitive
        coi_discrete_distribution[k+1] = beta_binom_density(k, n, alpha, beta)
    return coi_discrete_distribution

def beta_binom_density(k, n, alpha, beta):
    ''' return the probablity of occurence k for a beta binomial n, alpha, beta'''
    return 1.0*gamma(n+1)*gamma(alpha+k)*gamma(n+beta-k)*gamma(alpha+beta) / (gamma(k+1)*gamma(n-k+1)*gamma(alpha+beta+n)*gamma(alpha)*gamma(beta))

def beta_binom_mean(alpha,beta):
    mean = float(beta_distribution.stats(alpha,beta, moments='m'))
    return mean       

