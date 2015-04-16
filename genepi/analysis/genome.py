import os
import csv
import sys
import subprocess
import numpy as np
import pandas as pd

from genepi.genome import bp_per_cM

def plink_format(genomes):
    '''
    Transform dataframe into temporary PLINK formatted files (.map + .ped)

    By default, each line of the MAP file describes a single marker
    and must contain exactly 4 columns:
       chromosome (1-22, X, Y or 0 if unplaced)
       rs# or snp identifier
       Genetic distance (morgans)
       Base-pair position (bp units)

    '''

    def parse_SNP_id(id):
        '''
        Parse SNP ID from dataframe into .map formatted line:
        e.g. Pf.1.2402 --> 1 Pf.1.2402 0.0016013 2402
        '''
        _,chrom,pos=id.split('.')
        return "{0} {1} {2:.7f} {3}".format(chrom,id,float(pos)/bp_per_cM,pos)

    if not os.path.exists('output'):
        os.mkdir('output')
    with open(os.path.join('output','g.map'),'w') as f:
        for id in genomes.columns:
            line=parse_SNP_id(id)
            f.write(line+'\n')

    '''
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

    def parse_genome(id,genome):
        '''
        Parse genome ID and genome from dataframe into .ped formatted line:
        e.g. 0,00011... --> 0 g0 0 0 0 -9 1 1 1 1 1 1 2 2 2 2  ...
        '''
        return '0 g%d 0 0 0 -9 '%id + ' '.join([' '.join([str(x+1)]*2) for x in genome])

    with open(os.path.join('output','g.ped'),'w') as f:
        for id_genome in zip(genomes.index,genomes.values):
            line=parse_genome(*id_genome)
            f.write(line+'\n')

def ibd_finder():
    '''
    Call out to GERMLINE
    http://www.cs.columbia.edu/~gusev/germline/
    to find pairwise IBD segments
    '''
    with open('germline.stdin','r') as f:
        try:
            exe='./bin/germline'
            subprocess.check_call([exe,
                                   #'-haploid',
                                   #'-bin_out',
                                   '-min_m','10',
                                   '-bits','32',
                                   '-w_extend',
                                   '-err_hom','1',
                                   '-err_het','1'],
                                   stdin=f)
        except subprocess.CalledProcessError as e:
            codes={  1:'GERMLINE: complete',
                   -11:'GERMLINE: unhandled exception'}
            print('\n'+codes.get(e.returncode,'ERROR CODE=%d'%e.returncode))
        except OSError:
            print('GERMLINE executable not found at %s'%exe)

def cluster_finder():
    '''
    Call out to DASH
    http://www.cs.columbia.edu/~gusev/dash/
    to find clusters of sequences from IBD segments
    '''
    pass

def genome_analysis(file='simulations/GenomeReport.npz',reformat=True):
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

    cwd=os.path.dirname(os.path.realpath(__file__))
    has_file=lambda f: os.path.isfile(os.path.join(cwd,f))
    if reformat or not (has_file('g.map') and has_file('g.ped')):
        plink_format(genomes.iloc[-100:,:221]) # last hundred genomes, chrom-1

    ibd_finder()

    #cluster_finder()

if __name__ == '__main__':
    genome_analysis('../../examples/simulations/GenomeReport.npz',
                    reformat=True)
