import itertools
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

SNPs=defaultdict(list)

def get_N_SNPs():
    return sum([len(v) for v in SNPs.values()])

def from_txt_table(filename,allele_freqs=[]):
    '''
    Read SNP positions from file in following format:
    CHR    POS
    Pf3D7_01_v3    130339
    '''
    logger.info('Reading SNPs from file: %s' % filename)
    with open(filename) as f:
        for idx,content in enumerate(f.readlines()[1:]):
            chr,pos = content.split()
            freq=allele_freqs[idx] if allele_freqs else 0.5
            SNPs[int(chr.split('_')[1])].append((int(pos),freq))
    for c,s in SNPs.items():
        logger.debug('  Chr %d: SNPs (pos,freq): %s' % (c,s))

def get_recombination_locations(chrom):
    next_location=0
    locations=[]
    while next_location < P_falciparum_chromosomes_Mb[chrom]*bp_per_Mb:
        locations.append(next_location)
        d = int(math.ceil(random.expovariate(lambd=1.0/bp_per_morgan)))
        next_location+=d
    logger.debug('Chr %d: Recomb: %s' % (chrom,locations))
    return locations

class Genome:
    '''The discretized representation SNPs on chromosomes'''

    id=itertools.count()

    def __init__(self,genome=[]):
        self.id=Genome.id.next()
        if genome:
            self.genome=genome
        else:
            self.genome=random.getrandbits(get_N_SNPs())
        logger.debug('Genome: id=%d, genome=%s' % (self.id,self.print_bits()))

    def __repr__(self):
        return str(self.genome)

    def print_bits(self):
        return 'Genome: ' + ('{0:0>%db}' % get_N_SNPs()).format(self.genome)