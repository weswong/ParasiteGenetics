import random

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from genome import load_npz

def draw_weight_from_age(age):
    # new strains in symptomatic samples much higher weight
    # TODO: draw from EMOD DTK
    return random.uniform(0.2,1.0) if not age else random.uniform(0.01,0.2)

def sample_analysis(tx_file = 'simulations/TransmissionGeneticsReport.csv',
                    gn_file = 'simulations/GenomeReport.npz',
                    date_range = (15*365,20*365),
                    sample_rate = 0.01):
    try:
        tx=pd.read_csv(tx_file)
    except IOError as e:
        sys.exit(e)

    tx=tx.query('%d < day < %d'%date_range).set_index(['iid','day'])

    d=load_npz(gn_file)
    genomes=pd.DataFrame(d['genomes'],columns=d['header'])
    #print(genomes.head())

    samples={}
    for iid,group in tx.groupby(level=0):
        if random.random() < sample_rate:
            #print(group[['pid','gid']])
            if len(group) == 1:
                continue
            s = group.groupby(level=1).gid.apply(list)
            samples.update({iid:s.to_dict()})

    #print(samples)

    f,axs=plt.subplots(3,5,sharex=True,sharey=True,figsize=(15,8),num='Multi-strain infections')
    for i,(iid,gid_by_day) in enumerate(samples.items()):
        try:
            ax=axs[i//5,i%5]
        except:
            break
        print(gid_by_day)
        days = sorted(gid_by_day.keys())
        last_day = days[-1]
        weights=[]
        strains=[]
        for day in days:
            for gid in gid_by_day[day]:
                weights.append(draw_weight_from_age(last_day-day))
                strains.append(genomes.loc[gid].values)
        sum_weights = sum(weights)
        norm_weights = [w/sum_weights for w in weights]
        strains = np.array([s*w for (s,w) in zip(strains,norm_weights)])
        mean_alt = np.sum(strains, axis=0)
        logy=True
        ax.hist(mean_alt, bins=np.arange(0,1,0.01), log=logy, alpha=0.5, color='navy')
        #depths = np.random.randint(1,90,mean_alt.shape)
        depths = np.random.normal(80,20,mean_alt.shape)
        smear_alt = [np.random.binomial(max(5,d),max(0,min(1,m)))/float(d) for (m,d) in zip(mean_alt,depths)]
        ax.hist(smear_alt, bins=np.arange(0,1,0.01), log=logy, alpha=0.5, color='firebrick')
        ax.text(0.05,3e3,str(gid_by_day) + '\n' + str([ '%.2f' % elem for elem in norm_weights ]), fontsize=6)
        f.set_tight_layout(True)

if __name__ == '__main__':
    sample_analysis('../../examples/simulations/TransmissionGeneticsReport.csv',
                    '../../examples/simulations/GenomeReport.npz')
    plt.show()
