import matplotlib.pyplot as plt

import Genome as gn

gn.from_txt_table('barcode_loci.txt')

n_recombinations=[]
for _ in range(100):
    n=0
    for c in gn.P_falciparum_chromosomes_Mb.keys():
        if len(gn.get_recombination_locations(chrom=c))>1:
            n+=1 
    n_recombinations.append(n)
plt.hist(n_recombinations)
plt.show()