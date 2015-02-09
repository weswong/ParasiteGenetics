import itertools
import math
import random
from collections import defaultdict

import logging
logging.basicConfig(format='%(message)s')
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

bp_per_morgan = 1.5e6
bp_per_Mbp = 1e6
chrom_names=range(1,15) # + ['MT','Api']
chrom_lengths_Mbp=[0.643,0.947,1.1,1.2,1.3,1.4,1.4,1.3,1.5,1.7,2.0,2.3,2.7,3.3]
Pf_chrom_lengths=dict(zip(chrom_names,
                          [int(bp_per_Mbp*c) for c in chrom_lengths_Mbp]))
log.debug('Chromosome lengths:\n%s',Pf_chrom_lengths)

def initializeSNPs(SNP_source,bin_size=None):
    if SNP_source=='barcode':
        import barcode
        Genome.SNPs=barcode.SNPs
    else:
        raise Exception("Don't recognize SNP source type: %s", SNP_source)
    if bin_size:
        Genome.bin_size=bin_size
    else:
        min_dist=find_closest_SNPs(Genome.SNPs)
        Genome.bin_size=get_rounded_binning(min_dist,2)
    for snp in Genome.SNPs:
        snp.bin=int(snp.position/Genome.bin_size)
    log.debug(Genome.SNPs)

def find_closest_SNPs(snps):
    snp_pos_by_chrom=defaultdict(list)
    for snp in snps:
        snp_pos_by_chrom[snp.chromosome].append(snp.position)
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

def get_chromosome_bins(bp):
    return int(math.ceil(bp/Genome.bin_size))

def reference_genome():
    genome={}
    for chrom_name,chrom_len in Pf_chrom_lengths.iteritems():
        n_bins=get_chromosome_bins(chrom_len)
        chrom=[0]*n_bins
        genome[chrom_name]=chrom
    return genome

def iterate_SNPs():
    for s in Genome.SNPs:
        yield s.chromosome,s.bin,s.allele_freq

def get_n_SNPs():
    return len(Genome.SNPs)

def display_bit(b):
    #return str(b)
    return '*' if b else '-'

def pairwise(l):
    "l -> (l0,l1), (l2,l3), (l4, l5), ..."
    if len(l) % 2:
        l.append(None)
    a = iter(l)
    return itertools.izip(a, a)

def get_crossover_points(chrom,bp_per_morgan=bp_per_morgan):
    next_point=0
    xpoints=[]
    chrom_length=get_chromosome_bins(Pf_chrom_lengths[chrom])
    while next_point < chrom_length:
        if next_point:
            xpoints.append(next_point)
        d = int(math.ceil(random.expovariate(Genome.bin_size/bp_per_morgan)))
        next_point+=d
    #log.debug('Chr %s recomb: %s', chrom, xpoints)
    return xpoints

def crossover(c1,c2,xpoints):
    c3,c4=c1[:],c2[:]
    if not xpoints:
        return c3,c4
    for l1,l2 in pairwise(xpoints):
        c3[l1:l2], c4[l1:l2] = c4[l1:l2], c3[l1:l2]
    return c3,c4

def meiosis(in1,in2):
    #log.debug('Before meiosis:\n%s\n%s',in1,in2)
    genomes=[defaultdict(list) for _ in range(4)]
    for c in chrom_names:
        c1,c2=in1.genome[c],in2.genome[c]
        xpoints=get_crossover_points(chrom=c)
        #log.debug('Chr %d, xpoints=%s',c,xpoints)
        c3,c4=crossover(c1,c2,xpoints)
        outputs=sorted([c1,c2,c3,c4], key=lambda *args: random.random())
        for i in range(4):
            genomes[i][c]=outputs[i]
    out1,out2,out3,out4=(Genome(genomes[i]) for i in range(4))
    #log.debug('After meiosis:\n%s\n%s\n%s\n%s',out1,out2,out3,out4)
    return out1,out2,out3,out4

class SNP:
    '''
    The properties of a single nucleotide polymorphism
    '''

    def __init__(self,chrom,pos,freq=0.5,bin=None):
        self.chromosome=chrom
        self.position=pos
        self.allele_freq=freq
        self.bin=bin

    def __repr__(self):
        return 'SNP'+str((self.chromosome,self.position,self.bin))

class Genome:
    '''
    The discretized representation of SNPs on chromosomes
    '''

    id=itertools.count()
    SNPs=None
    bin_size=None

    def __init__(self,genome):
        self.id=Genome.id.next()
        log.debug('Genome: id=%d', self.id)
        self.genome=genome
        log.debug('%s', self)

    def __repr__(self):
        return self.display_barcode()
        #return self.display_genome()

    @classmethod
    def from_reference(cls):
        return cls(reference_genome())

    @classmethod
    def from_allele_frequencies(cls):
        genome=reference_genome()
        for c,b,f in iterate_SNPs():
            if random.random()<f:
                genome[c][b]=1
        return cls(genome)

    @classmethod
    def from_barcode(cls,barcode):
        genome=reference_genome()
        for (c,b,f),s in zip(iterate_SNPs(),barcode):
            genome[c][b]=s
        return cls(genome)

    def display_barcode(self):
        snp_values=[]
        for c,b,f in iterate_SNPs():
            snp_values.append(self.genome[c][b])
        return ''.join([display_bit(b) for b in snp_values])

    def display_genome(self):
        s=[]
        for chrom_name,chrom in self.genome.items():
            s.append('Chromosome %s' % chrom_name)
            s.append(''.join([display_bit(b) for b in chrom]))
        return '\n'.join(s)
