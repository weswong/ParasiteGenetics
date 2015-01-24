import os
import random
import numpy as np
import math
from itertools import izip, tee
import matplotlib.pyplot as plt
from plotting import adjust_spines
from scipy.stats import poisson

# TODO: smear number of variables sites based on independent draws per position from density-dependent detectability
# TODO: refine immune pressure function
# TODO: refine oocyst and hepatocyte multiplicity distributions
# TODO: refine gametocyte synchronicity dynamics

# Sandbox for mapping chokepoints in single parasite generation onto barcode diversity
def cotransmission(verbose=False):

    detected_variability=[]

    # Start with some diversity at all positions among sporozoites (S)
    S = [random.random() for k in range(24)]
    if verbose: print('S: %s' % S)

    for iteration in range(20):

        if verbose: print('ITERATION %d' % iteration)

        # Sample some infected hepatocytes (H) from sporozoites (S)
        H=[]
        for k in range( np.random.poisson(3) ):
            h = [ int(random.random() < b) for b in S ]
            H.append(h)
            if verbose: print('H%d: %s' % (k,h))
        if not H:
            if verbose: print('No infected hepatocytes.')
            break

        # Placeholder for some immunity-driven selection of certain strains
        # between rupturing hepatocytes (H) and high-density asexual stage (A)
        nH = len(H)
        fraction = 0.8
        A = H[:int(fraction*nH+0.5)]
        if verbose: print('unique A = %d' % len(A))
        if not A:
            if verbose: print('No viable blood-stage infection.')
            break

        # Some integrated asynchronicity between high-density asexual stages (A)
        # and gametocyte distributions (G) sampled by the mosquito
        G=np.zeros(24)
        norm=0
        for a in A:
            r = random.random()
            G += r*np.array(a)
            norm += r
        G /= norm
        if verbose: print('G: %s' % G)
        nvar = sum([0<b<1 for b in G])
        nvar_detected = sum([0.2<b<1 for b in G])
        if not nvar:
            if verbose: print('No more variable positions.')
            break
        else:
            detected_variability.append(nvar_detected)
            if verbose: print( 'Variable positions = %d, of which detectable %d' % (nvar,nvar_detected) )

        # Sampling new oocysts (O) from gametocyte distribution (G)
        O=[]
        for k in range( np.random.poisson(10) ):
            o = [ int(random.random() < b) for b in G ]
            O.append(o)
            if verbose: print('O%d: %s' % (k,o))
        if not O:
            if verbose: print('No oocyts.')
            break

        # Propagation of oocysts (O) to sporozoite distributions (S)
        S=np.zeros(24)
        for o in O:
            S += np.array(o)
        S /= len(O)
        if verbose: print('S: %s' % S)
        nvar = sum([0<b<1 for b in S])
        if not nvar:
            if verbose: print('No more variable positions.')
            break
        else:
            if verbose: print('Variable positions = %d' % nvar)

    if verbose: print('Generations of variability: %s' % detected_variability)
    return detected_variability

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

if __name__ == '__main__':

    variant_site_transitions = np.zeros((25,25), dtype='double')
    generations_of_variability = []

    for k in range(30000):
        v = cotransmission()
        #print(v)

        generations_of_variability.append(len(v))

        for i,j in pairwise(v):
            variant_site_transitions[i][j] += 1
        if len(v):
            variant_site_transitions[v[-1]][0] += 1 # for the final transition to zero

    #print('Generations of variability: %s' % generations_of_variability)
    #print('Variant-site number transitions: \n%s' % variant_site_transitions)

    fig=plt.figure('Generations of sustained variability', figsize=(8,6), facecolor='w')
    ax=fig.add_subplot(111)
    adjust_spines(ax,['left','bottom'])
    plt.hist(generations_of_variability, bins=np.arange(-0.5,10.5), color='dodgerblue')
    plt.xlim([-0.5,10.5])
    plt.xlabel('Generations of co-transmission to lose 24 variable sites')
    plt.tight_layout()
    #plt.savefig(os.path.join('..','figs','generations_of_sustained_variability.png'))

    def plot_transition(from_nvars):
        ax = plt.gca()
        dist = variant_site_transitions[from_nvars,:]
        dist /= sum(dist)
        plt.bar(np.arange(-0.5,24.5), dist, color='dodgerblue')
        adjust_spines(ax,['bottom'])
        plt.xlim([-0.5,12.5])
        ax.text(0.5,0.9,'Starting from %d variable sites' % from_nvars, transform=ax.transAxes, fontsize=12)
        k=np.arange(12)
        plt.plot(k, 0.18*(k==0)+0.82*poisson.pmf(k,from_nvars*0.45), color='firebrick', alpha=0.4)

    fig=plt.figure('Variant-site number transitions', figsize=(6,9), facecolor='w')
    ax=plt.subplot(611)
    plot_transition(10)

    plt.subplot(612)
    plot_transition(8)

    plt.subplot(613)
    plot_transition(6)

    plt.subplot(614)
    plot_transition(4)

    plt.subplot(615)
    plot_transition(2)

    plt.subplot(616)
    plot_transition(1)
    plt.gca().set_xlabel('# of variable sites in following generation')

    plt.tight_layout()
    #plt.savefig(os.path.join('..','figs','variant-site_number_transitions.png'))

    plt.show()