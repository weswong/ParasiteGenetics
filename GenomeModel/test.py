import random
from utils import log
import genome as gn

def add_reference(genomes):
    log.debug('\nReference:')
    genomes+=[gn.Genome.from_reference()]

def add_mutant(genomes):
    log.debug('\nComplete mutant from barcode:')
    genomes+=[gn.Genome.from_barcode([1]*gn.get_n_SNPs())]

def add_random(genomes,N):
  log.debug('\nRandom from allele frequencies:')
  genomes+=[gn.Genome.from_allele_frequencies() for _ in range(N)]

def init_test():
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    add_random(gg,N=10)

def meiosis_test(N):
    gg=[]
    add_reference(gg)
    add_mutant(gg)
    for _ in range(N):
       g1,g2 = random.sample(gg,2)
       gn.meiosis(g1,g2)


gn.initializeSNPs('barcode') ### why bin_size >= 2e5 crash?

init_test()

meiosis_test(5)
