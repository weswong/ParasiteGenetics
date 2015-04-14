import csv
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_unique=lambda x: len(x.unique())

def onward_transmissions(tx,ax1):
    n_transmit=tx.groupby('iidParent')['iid'].apply(n_unique).values
    infections=set(tx.iid.values)
    parents=set(tx.iidParent.values)
    zero_bin=len(infections^parents)
    n_transmit=np.append(n_transmit,[0]*zero_bin)
    ax1.hist(n_transmit,
             bins=np.arange(-0.5,10.5),normed=True,
             color='navy',alpha=0.2)
    ax1.set_xlabel('# onward transmissions')
    ax1.set_xlim([-1,11])

def genomes_per_infection(tx,ax2):
    n_genomes=tx.groupby('iid')['gid'].apply(n_unique).values
    ax2.hist(n_genomes,
             bins=np.arange(0.5,7.5),normed=True,
             color='navy',alpha=0.2)
    ax2.set_xlabel('# distinct genomes')
    ax2.set_xlim([0,8])

def repeat_genomes(tx):
    repeats=tx.groupby(['gid','day'])['pid'].count().unstack('day')
    repeats.dropna(thresh=3,inplace=True)
    f=plt.figure('RepeatGenomes')
    plt.imshow(repeats.values,interpolation='nearest',cmap='Reds')
    plt.xlabel('Timestep')
    plt.ylabel('Genome ID')
    f.set_tight_layout(True)

def transmission_analysis(file='simulations/TransmissionGeneticsReport.csv'):
    '''
    Analysis of the TransmissionGeneticsReport output
    '''
    try:
        tx=pd.read_csv(file)
    except IOError as e:
        sys.exit(e)
    f,(ax1,ax2)=plt.subplots(1,2,figsize=(10,5),num='InfectionProperties')
    onward_transmissions(tx,ax1)
    genomes_per_infection(tx,ax2)
    f.set_tight_layout(True)
    repeat_genomes(tx)

if __name__ == '__main__':
    transmission_analysis('../../examples/simulations/TransmissionGeneticsReport.csv')
    plt.show()
