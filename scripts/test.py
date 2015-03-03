import math
import random

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

from genepi import utils,report
import genepi.genome as gn
import genepi.infection as inf
import genepi.population as pop
import genepi.simulation as sim
import genepi.human as hn
import genepi.migration as mig

def add_reference(genomes):
    log.debug('REFERENCE')
    genomes+=[gn.Genome.from_reference()]

def add_mutant(genomes):
    log.debug('MUTANT')
    genomes+=[gn.Genome.from_barcode([1]*gn.get_n_SNPs())]

def add_random(genomes,N):
    log.debug('RANDOM')
    genomes+=[gn.Genome.from_allele_freq() for _ in range(N)]

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

def bitstring_test():
    A=[1 if random.random()<0.5 else 0 for _ in range(500)]
    s=''.join([str(b) for b in A])
    bs=int(s,2)
    bs2=sum(1<<i for i, b in enumerate(reversed(A)) if b)
    print(s,bs,bs2)

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
    print(utils.accumulate_cdf(A))

def transmit_test():
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
        #chosen=utils.choose_without_replacement(M,N)
        chosen=utils.choose_with_replacement(M,N)
        all_counts.update(chosen)
        log.debug(chosen)
    log.debug(all_counts)

def poisson_test(rate,n_tests):
    N=[]
    for _ in range(n_tests):
        N.append(utils.poissonRandom(rate))
    log.debug(N)

def population_test(tsteps):
    inf.log.setLevel(logging.DEBUG)
    pop.log.setLevel(logging.DEBUG)
    pop_params={'id':1,'parent':None,'n_humans':10,'n_infections':5}
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
    dt=sim.Params.sim_tstep
    for i in range(tsteps):
        G.send(dt)
        n=next(G)
        log.debug(n)

def reinfection_test(n_tsteps=5):
    sim.log.setLevel(logging.DEBUG)
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    i=inf.Infection(gg)
    t=0
    dt=21
    for _ in range(n_tsteps):
        txs=i.update(dt,0.1)
        for tx in txs:
            i.add_infection(tx)
        log.debug('t=%d\n%s'%(t,str(i)))
        t+=dt

def simulation_test():
    #sim.log.setLevel(logging.DEBUG)
    #pop.log.setLevel(logging.DEBUG)
    #inf.log.setLevel(logging.DEBUG)
    s=sim.Simulation()
    s.add_report(report.SimulationReport)
    #s.populations['Test'].add_infections([inf.Infection.from_random(1)])
    s.run()

def migration_test():
    sim.log.setLevel(logging.DEBUG)
    pop.log.setLevel(logging.DEBUG)
    inf.log.setLevel(logging.DEBUG)
    hn.log.setLevel(logging.DEBUG)
    mig.log.setLevel(logging.DEBUG)

    sim.Params.sim_duration = 21*6
    sim.Demographics.populations = {
        'Desert' : {
            'n_humans' : 10,
            'n_infections' : 0,
            'vectorial_capacity_fn' : sim.annual_cycle(0,coeff=0),
            'migration_rates' : {'Bog':0.1},
        },
        'Bog' : {
            'n_humans' : 10,
            'n_infections' : 2,
            'vectorial_capacity_fn' : sim.annual_cycle(0.1,coeff=0),
            'migration_rates' : {'Desert':0.1},
        }
    }
    s=sim.Simulation()
    s.run()

if __name__ == '__main__':
    gn.initializeSNPs('barcode')

    #SNP_test()
    #bitstring_test()
    #init_test()
    #meiosis_test(1)

    #plot_chokepoints()
    #accumulate_test()
    #poisson_test(0.9,100)
    #transmit_test()
    #reinfection_test()

    #sample_test(M=6,N=5,n_tests=10)
    #population_test(tsteps=2)
    #generator_test(tsteps=15)
    migration_test()

    #simulation_test()
