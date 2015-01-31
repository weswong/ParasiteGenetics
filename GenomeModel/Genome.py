import math
import random
from collections import defaultdict

from utils import logger

bp_per_morgan = 1.5e6
bp_per_Mb = 1e6

P_falciparum_chromosomes_Mb = { 1  : 0.643,
                                2  : 0.947,
                                3  : 1.1,
                                4  : 1.2,
                                5  : 1.3,
                                6  : 1.4,
                                7  : 1.4,
                                8  : 1.3,
                                9  : 1.5,
                                10 : 1.7,
                                11 : 2.0,
                                12 : 2.3,
                                13 : 2.7,
                                14 : 3.3,
                                #'MT' : ??,
                                #'Api': ??
                              }

SNP_positions=defaultdict(list)

def from_txt_table(filename):
    '''
    Read SNP positions from file in following format:
    CHR    POS
    Pf3D7_01_v3    130339
    '''
    logger.debug('Reading SNPs from file: %s' % filename)
    with open(filename) as f:
        for content in f.readlines()[1:]:
            chr,pos = content.split()
            SNP_positions[int(chr.split('_')[1])].append(int(pos))
    for c,p in SNP_positions.items():
        logger.debug('  Chr %d:\tSNPs: %s' % (c,p))

def get_recombination_locations(chrom):
    next_location=0
    locations=[]
    while next_location < P_falciparum_chromosomes_Mb[chrom]*bp_per_Mb:
        locations.append(next_location)
        d = int(math.ceil(random.expovariate(lambd=1.0/bp_per_morgan)))
        next_location+=d
    logger.debug('Chr %d:\tRecomb: %s' % (chrom,locations))
    return locations

class Genome:
    '''The discretized representation SNPs on chromosomes'''

    def __init__(self,genome):
        self.genome=genome

    @classmethod
    def random_from_allele_frequencies(cls):

        return cls(g)