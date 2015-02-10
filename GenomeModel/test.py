import random
import logging as log
import genome as gn
import infection as inf
import population as pop

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
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    add_random(gg,N=10)

def SNP_test():
    g=gn.reference_genome()
    for c,b,f in gn.iterate_SNPs():
        log.debug('Chrom %d (length %d): SNP at position %d' % (c,len(g[c]),b))

def meiosis_test(N):
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

def infection_test():
    inf.log.setLevel(log.DEBUG)
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    i=inf.Infection(gg)
    i.transmit()

def population_test():
    pop.log.setLevel(log.DEBUG)
    node={'id':1,'n_humans':100,'n_infections':10}
    p=pop.Population(**node)

gn.initializeSNPs('barcode')

#SNP_test()
#init_test()
#meiosis_test(1)

#plot_chokepoints()
#infection_test()
population_test()