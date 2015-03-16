import math
import random
import itertools
from collections import defaultdict
import hashlib

import logging
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import numpy as np # for fast meiosis operations on arrays

import utils

bp_per_morgan = 1.5e6
bp_per_Mbp    = 1e6

chrom_names       = range(1,15) # + ['MT','Api']
chrom_lengths_Mbp = [ 0.643, 0.947,  1.1,  1.2,   1.35,  1.42,  1.45,
                      1.5,   1.55,   1.7,  2.05,  2.3,   2.95,  3.3 ]
chrom_lengths_bp  = [ int(bp_per_Mbp*c) for c in chrom_lengths_Mbp ]
Pf_chrom_lengths  = dict(zip(chrom_names,chrom_lengths_bp))
log.debug('Chromosome lengths:\n%s',Pf_chrom_lengths)

def initialize_from(SNP_source,bin_size=None):

    SNPs=SNP.initialize_from(SNP_source)

    set_bin_size(SNPs,bin_size)

    set_chrom_breaks()

    set_binned_SNPs(SNPs)

def set_bin_size(SNPs,bin_size):
    if bin_size:
        Genome.bin_size_bp=bin_size
    else:
        min_dist=find_closest_SNPs(SNPs)
        Genome.bin_size_bp=get_rounded_binning(min_dist,2)

def find_closest_SNPs(SNPs):
    snp_pos_by_chrom=defaultdict(list)
    for snp in SNPs:
        snp_pos_by_chrom[snp.chrom].append(snp.pos)
    min_dist=[]
    for c,s in snp_pos_by_chrom.items():
        p=sorted(s)
        d=[j-i for i, j in itertools.izip(p[:-1], p[1:])]
        if d:
            min_dist=min(d+[min_dist])
    log.debug('Closest SNPs = %d bp', min_dist)
    return min_dist

def get_rounded_binning(min_dist, ndigits=1):
    if min_dist<1:
        raise Exception('Genome binning must be a postive integer.')
    log10_floor = math.floor(math.log10(min_dist))
    trunc_digits = -int(log10_floor) + (ndigits - 1)
    bin_size = round(min_dist - math.pow(10,-trunc_digits)//2, trunc_digits)
    log.debug('Rounded genome-discretization binning = %d bp', bin_size)
    return bin_size

def set_chrom_breaks():
    Genome.chrom_idxs={v:i for i,v in enumerate(chrom_names)}
    bins_per_chrom = [n_bins_from(bp) for bp in chrom_lengths_bp]
    Genome.chrom_breaks = list(utils.cumsum([0]+bins_per_chrom))
    #log.debug('chrom_breaks:',Genome.chrom_breaks)

def n_bins_from(bp):
    return int(math.ceil(bp/Genome.bin_size_bp))

def set_binned_SNPs(SNPs):
    Genome.SNP_bins=[]
    Genome.SNP_freqs=[]
    last_snp,last_bin,max_freq=(None,0,0)
    for snp in SNPs:
        snp_bin=int(snp.pos/Genome.bin_size_bp) + Genome.chrom_breaks[Genome.chrom_idxs[snp.chrom]]
        if last_snp and (snp.chrom,snp_bin)==(last_snp.chrom,last_bin):
            log.debug('overlap at chrom %d bin %d',snp.chrom,snp_bin)
            if snp.freq > max_freq:
                Genome.SNP_freqs[-1]=snp.freq # replace
                max_freq=snp.freq
        else:
            Genome.SNP_bins.append(snp_bin) # append
            Genome.SNP_freqs.append(snp.freq)
            last_snp,last_bin=(snp,snp_bin)

    log.info('%d of %d unique SNPs after discretization at %dbp binning',
             len(Genome.SNP_bins),len(SNPs),Genome.bin_size_bp)

def reference_genome():
    return np.zeros(genome_length(),dtype=np.uint8)

def num_SNPs():
    return len(Genome.SNP_bins)

def genome_length():
    return Genome.chrom_breaks[-1]

def display_bit(b):
    return '*' if b else '-'

def get_crossover_points(chrom_length,bp_per_morgan=bp_per_morgan):
    next_point=0
    xpoints=[]
    while next_point < chrom_length:
        if next_point:
            xpoints.append(next_point)
        d = int(math.ceil(random.expovariate(Genome.bin_size_bp/bp_per_morgan)))
        next_point+=d
    #log.debug('Chr %s recomb: %s', chrom, xpoints)
    return xpoints

# TODO: further optimization of crossover/meiosis as
#       size=(4,genome_length) matrix operations
#       faster option that only returns one byproduct?

# TODO: optimize distinct/hash check

# TODO: optimize barcode_as_long (i.e. faster serialize for reporting)

def crossover(c1,c2,xpoints):
    c3=np.copy(c1)
    c4=np.copy(c2)
    if not xpoints:
        return c3,c4
    for l1,l2 in utils.by_pair(xpoints):
        t = np.copy(c3[l1:l2])
        c3[l1:l2] = c4[l1:l2]
        c4[l1:l2] = t
    return c3,c4

#@profile
def meiosis(in1,in2):

    #log.debug('Before meiosis:\n%s\n%s',in1,in2)
    genomes=[reference_genome() for _ in range(4)]
    for idx,(start,end) in enumerate(utils.pairwise(Genome.chrom_breaks)):
        c1,c2=in1.genome[start:end],in2.genome[start:end]
        xpoints=get_crossover_points(chrom_length=len(c1))
        #log.debug('Chr %d, xpoints=%s',chrom_names[idx],xpoints)
        c3,c4=crossover(c1,c2,xpoints)
        outputs=sorted([c1,c2,c3,c4], key=lambda *args: random.random())
        for j in range(4):
            genomes[j][start:end]=outputs[j]
    out1,out2,out3,out4=(Genome(genomes[j]) for j in range(4))
    #log.debug('After meiosis:\n%s\n%s\n%s\n%s',out1,out2,out3,out4)
    return out1,out2,out3,out4

#@profile
def distinct(genomes):
    distinct=[]
    seen = set()
    for g in genomes:
        h=hash(g)
        if h not in seen:
            seen.add(h)
            distinct.append(g)
    return distinct

class SNP:
    '''
    The properties of a single nucleotide polymorphism
    '''

    def __init__(self,chrom,pos,freq=0.5,bin=None):
        self.chrom=chrom
        self.pos=pos
        self.freq=freq

    def __repr__(self):
        return 'genome.SNP(%d,%d)'%self.to_tuple()

    def to_tuple(self):
        return (self.chrom,self.pos)

    @staticmethod
    def initialize_from(SNP_source):
        if SNP_source=='barcode':
            from snp import barcode
            SNPs=barcode.SNPs
        elif SNP_source=='sequence':
            from snp import sequence
            SNPs=sequence.SNPs
        else:
            raise Exception("Don't recognize SNP source type: %s", SNP_source)
        return SNPs

class Genome:
    '''
    The discretized representation of SNPs on chromosomes
    '''

    id=itertools.count()

    chrom_idxs    = {} # chrom_break indices by chromosome name
    chrom_breaks  = [] # locations of chromosome breakpoints on genome
    SNP_bins      = [] # locations of variable positions on genome
    SNP_freqs     = [] # minor-allele frequency of binned SNPs
    bin_size_bp   = [] # base pairs per genome bin

    def __init__(self,genome):
        self.id=Genome.id.next()
        log.debug('Genome: id=%d', self.id)
        self.genome=genome
        log.debug('%s', self)

    def __repr__(self):
        return 'Genome(%s)'%self.genome

    def __str__(self):
        return self.display_barcode()
        #return self.display_genome()

    def __hash__(self):
        #return hash(tuple(self.genome))
        return int(hashlib.sha1(self.genome.view(np.uint8)).hexdigest(),16)

    @classmethod
    def from_reference(cls):
        return cls(reference_genome())

    @classmethod
    def from_allele_freq(cls):
        rands=np.random.random_sample((num_SNPs(),))
        barcode=rands<Genome.SNP_freqs
        return cls.from_barcode(barcode)

    @classmethod
    def from_barcode(cls,barcode):
        genome=reference_genome()
        np.put(genome,Genome.SNP_bins,barcode)
        return cls(genome)

    def barcode(self):
        return self.genome[Genome.SNP_bins]

    '''
    def barcode_as_long(self):
        return sum(1<<i for i,b in enumerate(reversed(self.barcode())) if b)
    '''

    def display_barcode(self):
        return ''.join([display_bit(b) for b in self.barcode()])

    def display_genome(self):
        s=[]
        for idx,(start,end) in enumerate(utils.pairwise(Genome.chrom_breaks)):
            s.append('Chromosome %s' % chrom_names[idx])
            s.append(''.join([display_bit(b) for b in self.genome[start:end]]))
        return '\n'.join(s)
