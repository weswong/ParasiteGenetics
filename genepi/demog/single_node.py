import logging
log = logging.getLogger(__name__)

import demog

def init():

    populations = {
        'Population #1' : {
            'n_humans' : 1350,
            'n_infections' : 500,
            'coi' : demog.init_coi_discrete_distribution(n=5, alpha=0.5, beta=8),
            'vectorial_capacity_fn': demog.calculate_r0(
                                                    {'r0' : [5,5,1.95,1.95, 2.35],'year' : [0,7*365,10*365,11*365,12*365],'regime': ['null'],'expire_rate' : 0.022},
                                                    #seasonal forcing term 
                                                    demog.annual_cycle(1.0,0.18,2), 
                                                    #average COI, + 1 because I shifted the distribution
                                                    demog.beta_binom_mean(alpha=0.5, beta=8) + 1
                                                    )
                            
        }
    }

    return populations
