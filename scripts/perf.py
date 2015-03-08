import random
import logging as log
import genepi.genome as gn

'''
Package:
  https://github.com/rkern/line_profiler

Installation:
  pip install line_profiler

Usage:
  Add @profile decorator to function
  kernprof -l -v perf.py

'''

@profile
def genome_speed_test():
    #gn.initializeSNPs('barcode')
    gn.initializeSNPs('sequence',bin_size=100)

    genomes=[gn.Genome.from_allele_freq() for _ in range(10)]

    for _ in range(1000):
        g1,g2 = random.sample(genomes,2)
        gn.meiosis(g1,g2)

genome_speed_test()
