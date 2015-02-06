from utils import log

import genome as gn
gn.Genome.initializeSNPs('barcode')

log.debug('\nRandom from allele frequencies:')
genomes=[gn.Genome.from_allele_frequencies() for _ in range(3)]
log.debug('\nReference:')
genomes+=[gn.Genome.from_reference()]
log.debug('\nComplete mutant from barcode:')
genomes+=[gn.Genome.from_barcode([1]*gn.Genome.get_n_SNPs())]

import numpy as np
import matplotlib.pyplot as plt
n_recombinations=[]
for _ in range(100):
    n=0
    for c in gn.Pf_chrom_lengths.index:
        if gn.get_recombination_locations(chrom=c):
            n+=1
    n_recombinations.append(n)
plt.hist(n_recombinations,bins=np.arange(-0.5,gn.Genome.get_n_SNPs()+0.5))
plt.show()

#import random
# def cxOnePoint(ind1, ind2):
#     size = min(len(ind1), len(ind2))
#     cxpoint = random.randint(1, size - 1)
#     ind1[cxpoint:], ind2[cxpoint:] = ind2[cxpoint:], ind1[cxpoint:]
#     return ind1, ind2
#
# g1=[1]*10
# g2=[0]*10
# cxOnePoint(g1,g2)
# print(g1,g2)
