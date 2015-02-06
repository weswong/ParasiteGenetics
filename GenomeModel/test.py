import genome as gn
gn.Genome.initializeSNPs('barcode')

print('\nRandom from allele frequencies:')
genomes=[gn.Genome.from_allele_frequencies() for _ in range(3)]
print('\nReference:')
genomes+=[gn.Genome.from_reference()]
print('\nComplete mutant from barcode:')
genomes+=[gn.Genome.from_barcode([1]*gn.Genome.get_n_SNPs())]

#import matplotlib.pyplot as plt
#n_recombinations=[]
#for _ in range(100):
#    n=0
#    for c in gn.P_falciparum_chromosomes_Mb.keys():
#        if len(gn.get_recombination_locations(chrom=c))>1:
#            n+=1
#    n_recombinations.append(n)
#plt.hist(n_recombinations)
#plt.show()

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
