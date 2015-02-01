import matplotlib.pyplot as plt

import Genome as gn

allele_freqs=[ 0.7567, 0.5334, 0.9282, 0.4096, 
               0.2809, 0.1045, 0.4659, 0.6340, 
               0.7400, 0.1839, 0.0155, 0.1021, 
               0.9507, 0.9131, 0.3911, 0.8800, 
               0.8641, 0.8476, 0.1940, 0.5053, 
               0.6744, 0.8951, 0.5995, 0.9217 ]
gn.from_txt_table('barcode_loci.txt',allele_freqs)

#n_recombinations=[]
#for _ in range(100):
#    n=0
#    for c in gn.P_falciparum_chromosomes_Mb.keys():
#        if len(gn.get_recombination_locations(chrom=c))>1:
#            n+=1 
#    n_recombinations.append(n)
#plt.hist(n_recombinations)
#plt.show()

genomes=[gn.Genome() for _ in range(10)]