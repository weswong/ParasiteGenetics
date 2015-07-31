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
        
def r0_data_prep(r0_params):
    for i in itertools.izip(r0_params['r0'], r0_params['r0'][1:]):
        if i[0] == i[1]:
            r0_params['regime'].append( 'constant' )
        else:
            r0_params['regime'].append('dynamic')
            
    return r0_params

def calculate_r0(r0_params, annual_cycle, mean_COI):
        
    r0_params = r0_data_prep(r0_params) 
    initial_r0 = r0_params['r0'][0]
    
    def r0_function(t, prob_mixed):
        index =bisect.bisect_left(r0_params['year'], t) 
        base_reproduction_rate = initial_r0 * r0_params['expire_rate'] * max(0, 1.0 - mean_COI/5.0) * annual_cycle(t)
        
        # if non-natural transmission rate is not constant'''
        def dynamic_r0(start_change_t, start_change_R0, end_change_t, end_change_R0, t):
            #end_r0 refers to the period in time where R0 stops changing and becomes constant again
            slope = (float(end_change_R0) - start_change_R0) / (end_change_t - start_change_t)
            delta_t = t-start_change_R0
            # if r0 is constant across the time points, should collapse down to r0 = base_reproduction_rate * r0_scale, where r0_scale = float(current_r0) / initial_r0

            r0 = base_reproduction_rate * (initial_r0 - (initial_r0 - start_change_R0) + (slope * delta_t)) / float(initial_r0)
            return r0
                
        if index == len(r0_params['year']):
            raise Exception('Cannot extrapolate beyond the given amount of time with which we have data for: maximum={max}, tried to do time={t}'.format(max=r0_params['year'][-1], t=t))
        else:
            r0 = dynamic_r0(r0_params['year'][index - 1], r0_params['r0'][index -1], 
                                r0_params['year'][index], r0_params['r0'][index],t)
        return r0
    return r0_function