import itertools
import math
import random
from utils import log
from collections import defaultdict

bp_per_morgan = 1.5e6
bp_per_Mbp = 1e6
chrom_names=range(1,15) # 'MT', 'Api'
chrom_lengths_Mbp=[0.643,0.947,1.1,1.2,1.3,1.4,1.4,1.3,1.5,1.7,2.0,2.3,2.7,3.3]
Pf_chrom_lengths=dict(zip(chrom_names,
                          [int(bp_per_Mbp*c) for c in chrom_lengths_Mbp]))
log.debug('Chromosome lengths:\n%s',Pf_chrom_lengths)

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

def display_bit(b):
    #return str(b)
    return '*' if b else '-'

@profile
def get_recombination_locations(chrom):
    next_location=0
    locations=[]
    chrom_length=Pf_chrom_lengths[chrom]
    while next_location < chrom_length:
        if next_location:
            locations.append(next_location)
        d = int(math.ceil(random.expovariate(lambd=1.0/bp_per_morgan)))
        next_location+=d
    #log.debug('Chr %s recomb: %s', chrom, locations)
    return locations

class SNP:
    '''
    An object containing properties of a single nucleotide polymorphism
    '''
    def __init__(self,chrom,pos,freq=0.5,bin=None):
        self.chromosome=chrom
        self.position=pos
        self.allele_freq=freq
        self.bin=bin

    def __repr__(self):
        return 'SNP'+str((self.chromosome,self.position,self.bin))

class Genome:
    '''The discretized representation of SNPs on chromosomes'''

    id=itertools.count()
    SNPs=None
    bin_size=None

    @classmethod
    @profile
    def initializeSNPs(cls,SNP_source,bin_size=[]):
        if SNP_source=='barcode':
            import barcode
            cls.SNPs=barcode.SNPs
        else:
            raise Exception("Don't recognize SNP source type: %s", SNP_source)
        if bin_size:
            cls.bin_size=bin_size
        else:
            min_dist=find_closest_SNPs(cls.SNPs)
            cls.bin_size=get_rounded_binning(min_dist,2)
        for snp in cls.SNPs:
            snp.bin=int(snp.position//cls.bin_size)
        log.debug(cls.SNPs)

    def __init__(self,genome):
        self.id=Genome.id.next()
        log.debug('Genome: id=%d', self.id)
        self.genome=genome
        log.debug('%s', self)

    def __str__(self):
        return self.display_barcode()
        #return self.display_genome()

    @classmethod
    def reference_genome(cls):
        genome={}
        for chrom_name,chrom_len in Pf_chrom_lengths.iteritems():
            n_bins=int(chrom_len/cls.bin_size)
            chrom=[0]*n_bins
            genome[chrom_name]=chrom
        return genome

    @classmethod
    def from_reference(cls):
        return cls(cls.reference_genome())

    @classmethod
    @profile
    def from_allele_frequencies(cls):
        genome=cls.reference_genome()
        for c,b,f in cls.iterate_SNPs():
            if random.random()<f:
                genome[c][b]=1
        return cls(genome)

    @classmethod
    def from_barcode(cls,barcode):
        genome=cls.reference_genome()
        for (c,b,f),s in zip(cls.iterate_SNPs(),barcode):
            genome[c][b]=s
        return cls(genome)

    @classmethod
    def iterate_SNPs(cls):
        for s in cls.SNPs:
            yield s.chromosome,s.bin,s.allele_freq

    def display_barcode(self):
        snp_values=[]
        for c,b,f in self.iterate_SNPs():
            snp_values.append(self.genome[c][b])
        return ''.join([display_bit(b) for b in snp_values])

    def display_genome(self):
        s=[]
        for chrom_name,chrom in self.genome.items():
            s.append('Chromosome %s' % chrom_name)
            s.append(''.join([display_bit(b) for b in chrom]))
        return '\n'.join(s)

    @classmethod
    def get_n_SNPs(cls):
        return len(cls.SNPs)
