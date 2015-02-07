from utils import log

import genome as gn
gn.Genome.initializeSNPs('barcode')

log.debug('\nRandom from allele frequencies:')
genomes=[gn.Genome.from_allele_frequencies() for _ in range(3)]
log.debug('\nReference:')
genomes+=[gn.Genome.from_reference()]
log.debug('\nComplete mutant from barcode:')
genomes+=[gn.Genome.from_barcode([1]*gn.Genome.get_n_SNPs())]

for c in gn.Pf_chrom_lengths.keys():
    r=gn.get_recombination_locations(chrom=c)

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
