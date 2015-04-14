import sys
import subprocess
import numpy as np
import pandas as pd

def plink_format(genomes):
    '''
    Transform dataframe into temporary PLINK formatted files (.map + .ped)

    By default, each line of the MAP file describes a single marker
    and must contain exactly 4 columns:
       chromosome (1-22, X, Y or 0 if unplaced)
       rs# or snp identifier
       Genetic distance (morgans)
       Base-pair position (bp units)

    The PED file is a white-space (space or tab) delimited file;
    the first six columns are mandatory:
       Family ID
       Individual ID
       Paternal ID
       Maternal ID
       Sex (1=male; 2=female; other=unknown)
       Phenotype

    Genotypes (column 7 onwards) should also be white-space delimited;
    they can be any character (e.g. 1,2,3,4 or A,C,G,T or anything else)
    except 0 which is, by default, the missing genotype character.
    All markers should be biallelic. All SNPs (whether haploid or not)
    must have two alleles specified. Either Both alleles should be
    missing (i.e. 0) or neither.

    '''
    print(genomes[:10])

def ibd_finder():
    '''
    Call out to GERMLINE
    http://www.cs.columbia.edu/~gusev/germline/
    to find pairwise IBD segments
    '''
    pass

def cluster_finder():
    '''
    Call out to DASH
    http://www.cs.columbia.edu/~gusev/dash/
    to find clusters of sequences from IBD segments
    '''
    pass

def genome_analysis(file='simulations/GenomeReport.npz'):
    '''
    Analysis of the GenomeReport output
    '''
    try:
        with np.load(file) as data:
            A = data['genomes']
            header=data['header']
    except IOError as e:
        sys.exit(e)

    genomes=pd.DataFrame(A,columns=header)
    plink_format(genomes)
    ibd_finder()
    cluster_finder()

if __name__ == '__main__':
    genome_analysis('../../examples/simulations/GenomeReport.npz')
