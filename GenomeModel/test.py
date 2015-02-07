import random
from utils import log
import genome as gn
gn.initializeSNPs('barcode') ### why bin_size >= 2e5 crash?

genomes=[]

#log.debug('\nRandom from allele frequencies:')
#genomes+=[gn.Genome.from_allele_frequencies() for _ in range(3)]

log.debug('\nReference:')
genomes+=[gn.Genome.from_reference()]

log.debug('\nComplete mutant from barcode:')
genomes+=[gn.Genome.from_barcode([1]*gn.get_n_SNPs())]

for _ in range(5):
   g1,g2 = random.sample(genomes,2)
   gn.meiosis(g1,g2)
