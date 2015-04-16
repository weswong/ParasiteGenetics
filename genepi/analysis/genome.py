import os
import csv
import sys
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from genepi.genome import bp_per_cM

cwd=os.path.dirname(os.path.realpath(__file__))

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

    file=os.path.join('output','plink.map')
    print('Writing file: %s'%file)
    with open(file,'w') as f:
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

    def parse_genome(id,genome=[]):
        '''
        Parse genome ID and genome from dataframe into .ped formatted line:
        e.g. 0,00011... --> 0 g0 0 0 0 -9 1 1 1 1 1 1 2 2 2 2  ...
        '''
        return '0 g%d 0 0 0 -9 '%id + ' '.join([' '.join([str(x+1)]*2) for x in genome])

    file=os.path.join('output','plink.ped')
    print('Writing file: %s'%file)
    with open(file,'w') as f:
        for id_genome in zip(genomes.index,genomes.values):
            line=parse_genome(*id_genome)
            f.write(line+'\n')

    file=os.path.join('output','plink.fam')
    print('Writing file: %s'%file)
    with open(file,'w') as f:
        for id in genomes.index:
            line=parse_genome(id)
            f.write(line+'\n')

def ibd_finder():
    '''
    Call out to GERMLINE
    http://www.cs.columbia.edu/~gusev/germline/
    to find pairwise IBD segments
    '''
    with open(os.path.join(cwd,'germline.stdin'),'r') as f:
        try:
            exe=os.path.join(cwd,'bin','germline')
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
            codes={1:'GERMLINE: complete'}
            print(codes.get(e.returncode,'\nERROR CODE=%d'%e.returncode))
        except OSError:
            print('GERMLINE executable not found at %s'%exe)

def cluster_finder():
    '''
    Call out to DASH
    http://www.cs.columbia.edu/~gusev/dash/
    to find clusters of sequences from IBD segments
    '''
    exe=os.path.join(cwd,'bin','dash_cc')
    cmds=['cat output/germline.match',
          'cut -f 1,2,4',
          '%s output/plink.fam output/dash'%exe]
    try:
        output=subprocess.check_output('|'.join(cmds), shell=True)
    except subprocess.CalledProcessError as e:
        codes={1:'DASH: complete'}
        print(codes.get(e.returncode,'\nERROR CODE=%d'%e.returncode))
    except OSError:
        print('DASH executable not found at %s'%exe)

def ibd_analysis():
    df=pd.read_csv('output/germline.match',
                   delim_whitespace=True,
                   names=['famId1','indId1','famId2','indId2',
                          'chrom','start','end','startSNP','endSNP',
                          'bits','dist','unit','mismatches','homo1','homo2'])

    # IBD segment lengths
    f=plt.figure('SegmentLengthIBD')
    df.dist.plot(kind='hist',bins=30,color='navy',alpha=0.2,fig=f)
    plt.xlabel('IBD segment length (cM)')

    # IBD segment locations + frequency
    counts=df.groupby(['chrom','start','end'])['dist'].count()
    nchrom=len(counts.index.levels[0])
    f,axs=plt.subplots(nchrom,1,num='SegmentMapIBD')
    for ichrom,(chrom,ibd_counts) in enumerate(counts.groupby(level=0)):
        ax=axs[ichrom] if nchrom>1 else axs
        for idx,(rng,cnt) in enumerate(ibd_counts[chrom].iteritems()):
            ax.plot([x/1e6 for x in rng],[idx]*2,linewidth=0.2+1.2*np.log10(cnt),c='navy',alpha=0.8)
            ax.set_title('Chromosome %s (MB)'%chrom,x=0.12,y=0.85,fontsize=12)
            ax.get_yaxis().set_visible(False)
    f.set_tight_layout(True)

    # IBD fractions
    shared=df.groupby(['indId1','indId2'])['dist'].sum()
    f=plt.figure('SharedLengthIBD')
    shared.plot(kind='hist',bins=30,color='navy',alpha=0.2,fig=f)
    plt.xlabel('Total IBD length (cM)')

    # pairwise chromosome painter
    g1,g2=shared.argmax()
    print('Most similar genomes: %s and %s (IBD=%0.1f cM)'%(g1,g2,shared[(g1,g2)]))

    # IBD fraction network

def cluster_analysis():
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
    has_file=lambda f: os.path.isfile(os.path.join(cwd,'output',f))

    if reformat or not (has_file('plink.map') and has_file('plink.ped')):
        # chrom-1: 221
        # chrom-2: 593
        # chrom-3: 1023
        # chrom-4: 1535
        # TODO: funny behavior with GERMLINE running on multiple chromosomes
        plink_format(genomes.iloc[::10,1023:1535]) # 10x downsampling, chrom-3/4

    if reformat or not (has_file('germline.match') and has_file('plink.fam')):
        ibd_finder()
    ibd_analysis()

    # TODO: handling of multiple chromosomes?
    # if reformat or not has_file('dash.clst'):
    #     cluster_finder()
    # cluster_analysis()

if __name__ == '__main__':
    genome_analysis('../../examples/simulations/GenomeReport.npz',
                    reformat=False)
    plt.show()
