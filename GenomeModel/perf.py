from utils import log
import genome as gn

@profile
def genome_speed_test():
    gn.Genome.initializeSNPs('barcode')

    genomes=[gn.Genome.from_allele_frequencies() for _ in range(10000)]

    for g in genomes:
        for c in gn.Pf_chrom_lengths.keys():
            r=gn.get_recombination_locations(chrom=c)

#@profile
def test_dict_lookup(c):
    L=dict(zip(gn.chrom_names,gn.chrom_lengths_Mbp))
    l=L[c]
    print(l)

#@profile
def test_pandas_lookup(c):
    L=gn.Pf_chrom_lengths
    l=L.loc[c]
    print(l)

def build_SNP_class_list():
    snps=[]
    gn.Genome.initializeSNPs('barcode')
    for c,p,b,f in gn.Genome.SNPs[['chromosome', 'position', 'bin', 'allele_freq']].values:
        snps.append(gn.SNP(c,p,f,b))
    return snps

#@profile
def test_pandas_iteration():
    for c,b,f in gn.Genome.iterate_SNPs():
        pass

#@profile
def test_SNP_class_iteration(snps):
    for s in snps:
        c,b,f=s.chromosome,s.bin,s.allele_freq

def pandas_vs_dict():
    import random
    rands=[random.choice(gn.chrom_names) for _ in range(100)]
    for r in rands:
        test_pure_lookup(r)
        test_pandas_lookup(r)

def pandas_DF_vs_SNP_class():
    snps=build_SNP_class_list()
    test_pandas_iteration()
    test_SNP_class_iteration(snps)

#pandas_vs_dict()
genome_speed_test()
#pandas_DF_vs_SNP_class()
