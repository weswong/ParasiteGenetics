import random
import genome as gn
import infection as inf
import population as pop
import simulation as sim

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def add_reference(genomes):
    log.debug('REFERENCE')
    genomes+=[gn.Genome.from_reference()]

def add_mutant(genomes):
    log.debug('MUTANT')
    genomes+=[gn.Genome.from_barcode([1]*gn.get_n_SNPs())]

def add_random(genomes,N):
    log.debug('RANDOM')
    genomes+=[gn.Genome.from_allele_frequencies() for _ in range(N)]

def init_test():
    gn.log.setLevel(logging.DEBUG)
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    add_random(gg,N=10)

def SNP_test():
    g=gn.reference_genome()
    for c,b,f in gn.iterate_SNPs():
        log.debug('Chrom %d (length %d): SNP at position %d' % (c,len(g[c]),b))

def meiosis_test(N):
    gn.log.setLevel(logging.DEBUG)
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    for _ in range(N):
        log.debug('MEIOSIS')
        g1,g2 = random.sample(gg,2)
        gn.meiosis(g1,g2)

def plot_chokepoints():
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns

    o=[inf.sample_n_oocysts() for _ in range(2000)]
    h=[inf.sample_n_hepatocytes() for _ in range(2000)]
    data=pd.DataFrame({'oocysts':o,'hepatocytes':h})

    sns.set(style='ticks')

    f,(ax1,ax2)=plt.subplots(2,1,sharex=True)
    ax1.hist(o,bins=np.arange(-0.5,49.5),alpha=0.3)
    ax1.set_title('Oocysts')
    ax2.hist(h,bins=np.arange(-0.5,49.5),alpha=0.3)
    ax2.set_title('Hepatocytes')
    ax2.set_xlim([-1,40])

    sns.jointplot(x='oocysts',y='hepatocytes',data=data,
                  kind='scatter',stat_func=None,dropna=False,color='navy',
                  marginal_kws={'bins':np.arange(-0.5,49.5)},
                  joint_kws={'alpha':0.1},
                  xlim=[-1,40],ylim=[-1,40])

    plt.show()

def accumulate_test():
    #A=[random.random() for _ in range(20)]
    A=[1]*random.randint(1,10)
    print(inf.accumulate_cdf(A))

def infection_test():
    inf.log.setLevel(logging.DEBUG)
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    i=inf.Infection(gg)
    i.transmit()

def sample_test(M,N,n_tests):
    from collections import Counter
    all_counts=Counter()
    for _ in range(n_tests):
        log.debug('Choosing %d from %d',M,N)
        chosen=pop.choose_without_replacement(M,N)
        all_counts.update(chosen)
        log.debug(chosen)
    log.debug(all_counts)

def population_test(tsteps):
    inf.log.setLevel(logging.DEBUG)
    pop.log.setLevel(logging.DEBUG)
    pop_params={'id':1,'n_humans':10,'n_infections':5}
    p=pop.Population(**pop_params)
    for tstep in range(tsteps):
        transmitted_infections=[]
        for i in p.infections.values():
            if random.random() < 0.5:
                transmitted_infections.append(i.transmit())
        p.add_infections(transmitted_infections)
        log.debug('  %s',p)

def generator_test(tsteps):
    log.debug('GENERATOR')
    G=sim.Params.infectious_generator()
    G.send(None)
    for i in range(tsteps):
        G.send(sim.Params.sim_tstep)
        n=next(G)
        log.debug(n)

def simulation_test():
    sim.log.setLevel(logging.DEBUG)
    pop.log.setLevel(logging.DEBUG)
    #inf.log.setLevel(logging.DEBUG)
    s=sim.Simulation()
    #s.populations['Test'].add_infections([inf.Infection.from_random(1)])
    s.run()

def migration_test():
    sim.log.setLevel(logging.DEBUG)
    pop.log.setLevel(logging.DEBUG)
    #inf.log.setLevel(logging.DEBUG)
    sim.Params.sim_duration = 63
    sim.Demographics.populations = {
        'Test1' : {
            'n_humans' : 1,
            'n_infections' : 1,
            'vectorial_capacity' : sim.annual_cycle(0,coeff=0),
            'migration_rates' : {'Test2':0.1},
        },
        'Test2' : {
            'n_humans' : 1,
            'n_infections' : 0,
            'vectorial_capacity' : sim.annual_cycle(0.3,coeff=0),
            'migration_rates' : {'Test1':0},
        }
    }
    s=sim.Simulation()
    s.run()

if __name__ == '__main__':
    gn.initializeSNPs('barcode')

    #SNP_test()
    #init_test()
    #meiosis_test(1)

    #plot_chokepoints()
    #accumulate_test()
    #infection_test()

    sample_test(M=2,N=5,n_tests=1000)
    #population_test(tsteps=2)
    #generator_test(tsteps=5)
    #simulation_test()
    #migration_test()
