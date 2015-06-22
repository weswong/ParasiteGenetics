import math
import bisect
import random
import itertools
from collections import defaultdict
import hashlib
import sys

import logging
log = logging.getLogger(__name__)

import numpy as np # for fast meiosis operations on arrays

import utils
from snp.snp import SNP
from transmission import Transmission

bp_per_morgan = 1.5e6
bp_per_cM     = bp_per_morgan/1e2
bp_per_Mbp    = 1e6

# mutation_prob = 9.7e-9 / 2.0       #Bopp et al, 2013. Miottoic evolution of P falciparum shows a stable core genome but recombination in antigen families  per day, per position

chrom_names       = range(1,15) # + ['MT','Api']
chrom_lengths_Mbp = [ 0.643, 0.947,  1.1,  1.2,   1.35,  1.42,  1.45,
                      1.5,   1.55,   1.7,  2.05,  2.3,   2.95,  3.3 ]
chrom_lengths_bp  = [ int(bp_per_Mbp*c) for c in chrom_lengths_Mbp ]
Pf_chrom_lengths  = dict(zip(chrom_names,chrom_lengths_bp))
log.debug('Chromosome lengths:\n%s',Pf_chrom_lengths)

def initialize_from(SNP_source,bin_size=None,min_allele_freq=0):
    SNPs=SNP.initialize_from(SNP_source,min_allele_freq)
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

def snp_bin_from_chrom_pos(chrom,pos):
    return int(pos/Genome.bin_size_bp) + Genome.chrom_breaks[Genome.chrom_idxs[chrom]]

def set_binned_SNPs(SNPs):
    Genome.SNP_bins=[]
    Genome.SNP_freqs=[]
    Genome.SNP_names=[]
    last_snp,last_bin,max_freq=(None,0,0)
    for snp in SNPs:
        snp_bin=snp_bin_from_chrom_pos(snp.chrom,snp.pos)
        if last_snp and (snp.chrom,snp_bin)==(last_snp.chrom,last_bin):
            log.debug('overlap at chrom %d bin %d',snp.chrom,snp_bin)
            if snp.freq > max_freq:
                Genome.SNP_freqs[-1]=snp.freq # replace
                Genome.SNP_names[-1]='Pf.%d.%d'%(snp.chrom,snp.pos)
                max_freq=snp.freq
        else:
            Genome.SNP_bins.append(snp_bin) # append
            Genome.SNP_freqs.append(snp.freq)
            Genome.SNP_names.append('Pf.%s.%s'%(snp.chrom,snp.pos))
            last_snp,last_bin=(snp,snp_bin)

    Genome.bin_fitness = np.ones(genome_length())

    log.info('%d of %d unique SNPs after discretization at %dbp binning',
             len(Genome.SNP_bins),len(SNPs),Genome.bin_size_bp)

def add_locus(chrom,pos,name,fitness,freq=0):
    add_bin = snp_bin_from_chrom_pos(chrom,pos)
    if add_bin in Genome.SNP_bins:
        idx = bisect.bisect_left(Genome.SNP_bins, add_bin)
        log.warning('SNP already exists at bin %d: chrom,pos=(%d,%d)' % (add_bin, chrom, pos))
        Genome.SNP_names[idx] = name
        Genome.SNP_freqs[idx] = freq
    else:
        Genome.SNP_bins.append(add_bin)
        Genome.SNP_names.append(name)
        Genome.SNP_freqs.append(freq)
    Genome.bin_fitness[add_bin] = fitness

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

def meiosis(in1,in2,N=4):
    if N > 4:
        raise IndexError('Maximum of four distinct meiotic products to sample.')
    genomes=[reference_genome() for _ in range(4)]
    for idx,(start,end) in enumerate(utils.pairwise(Genome.chrom_breaks)):
        c1,c2=in1.genome[start:end],in2.genome[start:end]
        xpoints=get_crossover_points(chrom_length=len(c1))
        #log.debug('Chr %d, xpoints=%s',chrom_names[idx],xpoints)
        c3,c4=crossover(c1,c2,xpoints)
        outputs=sorted([c1,c2,c3,c4], key=lambda *args: random.random())
        for j in range(4):
            genomes[j][start:end]=outputs[j]
    return [Genome(genomes[j]) for j in range(N)]

def single_meiotic_product(in1,in2):
    out_genome=reference_genome()
    for idx,(start,end) in enumerate(utils.pairwise(Genome.chrom_breaks)):
        R=random.random()
        if R<0.25:
            out_genome[start:end]=in1.genome[start:end]
        elif R<0.5:
            out_genome[start:end]=in2.genome[start:end]
        else:
            xpoints=get_crossover_points(chrom_length=(end-start))
            #log.debug('Chr %d, xpoints=%s',chrom_names[idx],xpoints)
            c3,c4=crossover(in1.genome[start:end],in2.genome[start:end],xpoints)
            out_genome[start:end]=c3 if R<0.75 else c4
    return Genome(out_genome)

def distinct_sporozoites_from(gametocyte_pairs,n_products):
    transmitted_sporozoites=[]
    for (g1,g2),N in zip(gametocyte_pairs,n_products):
        if g1.id==g2.id:
            #log.debug('Selfing of gametocytes (id=%d)\n%s',g1.id,g1)
            t=Transmission((g1.id,g2.id),g1)
            transmitted_sporozoites.append(t)
            continue
        if N>1:
            meiotic_products=meiosis(g1,g2,N)
            #log.debug('Meiosis: %s',[str(mp) for mp in meiotic_products])
            tt=[Transmission((g1.id,g2.id),g) for g in meiotic_products]
            transmitted_sporozoites.extend(tt)
        else:
            selected_product=single_meiotic_product(g1,g2)
            #log.debug('Single meiotic product: %s',selected_product)
            t=Transmission((g1.id,g2.id),selected_product)
            transmitted_sporozoites.append(t)
    return distinct(transmitted_sporozoites,id_fn=lambda t:t.genome.id)


def distinct(genomes,
             id_fn=lambda g:g.id):
             #id_fn=lambda g: hash(g)):
    distinct=[]
    seen = set()
    for g in genomes:
        h=id_fn(g)
        if h not in seen:
            seen.add(h)
            distinct.append(g)
    return distinct

class Genome:
    '''
    The discretized representation of SNPs on chromosomes
    '''

    id=itertools.count()
    hash_to_id   = {}

    chrom_idxs    = {} # chrom_break indices by chromosome name
    chrom_breaks  = [] # locations of chromosome breakpoints on genome
    SNP_bins      = [] # locations of variable positions on genome
    SNP_names     = [] # chrom.pos encoding of binned SNPs
    SNP_freqs     = [] # minor-allele frequency of binned SNPs
    bin_fitness   = [] # relative fitness at each binned site
    bin_size_bp   = [] # base pairs per genome bin

    # TODO: find a better thread-safe way of letting Genome know what
    #       Simulation to notify on reportable events
    sim = None
    @classmethod
    def set_simulation_ref(cls,sim):
        cls.sim=sim

    def __init__(self,genome,mod_fns=[]):
        self.genome=genome
        for fn in mod_fns:
            fn(self.genome)
        h=hash(self)
        id=Genome.hash_to_id.get(h,None)
        #id=None
        if id:
            self.id=id
        else:
            self.id=Genome.id.next()
            Genome.hash_to_id[h]=self.id
            if Genome.sim:
                Genome.sim.notify('genome.init',self)
        log.debug('Genome: id=%d', self.id)
        log.debug('%s', self)

    def __repr__(self):
        return 'Genome(%s)'%self.genome

    def __str__(self):
        return self.display_barcode()
        #return self.display_genome()

    def __hash__(self):
        h=hashlib.sha1
        #h=hashlib.md5
        return int(h(self.genome.view(np.uint8)).hexdigest(),16)

    @classmethod
    def from_reference(cls):
        return cls(reference_genome())

    @classmethod
    def from_allele_freq(cls,mod_fns=[]):
        rands=np.random.random_sample((num_SNPs(),))
        barcode=rands<Genome.SNP_freqs
        return cls.from_barcode(barcode,mod_fns)

    @classmethod
    def from_barcode(cls,barcode,mod_fns=[]):
        genome=reference_genome()
        np.put(genome,Genome.SNP_bins,barcode)
        return cls(genome,mod_fns)

    def fitness(self):
        m=self.bin_fitness[self.genome!=0] # NB: assuming binary SNPs
        return np.product(m) if m.size else 1.

    def barcode(self,sites=None):
        if not sites:
            sites=Genome.SNP_bins
        return self.genome[sites]

    def display_barcode(self):
        return ''.join([display_bit(b) for b in self.barcode()])

    def display_genome(self):
        s=[]
        for idx,(start,end) in enumerate(utils.pairwise(Genome.chrom_breaks)):
            s.append('Chromosome %s' % chrom_names[idx])
            s.append(''.join([display_bit(b) for b in self.genome[start:end]]))
        return '\n'.join(s)
