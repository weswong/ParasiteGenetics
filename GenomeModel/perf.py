import random
from utils import log
import genome as gn

@profile
def genome_speed_test():
    gn.initializeSNPs('barcode')

    genomes=[gn.Genome.from_allele_frequencies() for _ in range(100)]

    for _ in range(100):
        g1,g2 = random.sample(genomes,2)
        gn.meiosis(g1,g2)

genome_speed_test()
