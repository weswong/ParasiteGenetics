import os
import csv
import glob
import sys
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch as fpatch
import networkx as nx

try:
    import seaborn as sns
    import matplotlib.colors
    sns.set_style("white",{'grid.color': '.93'})
    sns.set_context("notebook")
    cmap = matplotlib.colors.ListedColormap(sns.color_palette("husl",10))
except:
    print('Install seaborn package for more pleasing aesthetics.')
    cmap='hsv'

from genepi.genome import bp_per_cM,chrom_names,chrom_lengths_Mbp

cwd=os.path.dirname(os.path.realpath(__file__))

def plink_format(genomes,chrom_name=''):
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

    file=os.path.join('output','plink%s.map'%chrom_name)
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

    file=os.path.join('output','plink%s.ped'%chrom_name)
    print('Writing file: %s'%file)
    with open(file,'w') as f:
        for id_genome in zip(genomes.index,genomes.values):
            line=parse_genome(*id_genome)
            f.write(line+'\n')

    file=os.path.join('output','plink%s.fam'%chrom_name)
    print('Writing file: %s'%file)
    with open(file,'w') as f:
        for id in genomes.index:
            line=parse_genome(id)
            f.write(line+'\n')

def ibd_finder(chrom_name=''):
    '''
    Call out to GERMLINE
    http://www.cs.columbia.edu/~gusev/germline/
    to find pairwise IBD segments
    '''
    with open(os.path.join(cwd,'germline.stdin'),'w') as f:
        f.write('1\noutput/plink%s.map\noutput/plink%s.ped\noutput/germline%s'%(chrom_name,chrom_name,chrom_name))
    with open(os.path.join(cwd,'germline.stdin'),'r') as f:
        try:
            exe=os.path.join(cwd,'bin','germline')
            subprocess.check_call([exe,
                                   #'-haploid',
                                   #'-bin_out',
                                   '-min_m','10',
                                   '-bits','32',
                                   '-w_extend',
                                   '-err_hom','0',
                                   '-err_het','0'],
                                   stdin=f)
        except subprocess.CalledProcessError as e:
            codes={1:'GERMLINE: complete'}
            print(codes.get(e.returncode,'\nERROR CODE=%d'%e.returncode))
        except OSError:
            print('GERMLINE executable not found at %s'%exe)

def cluster_finder(chrom_name=''):
    '''
    Call out to DASH
    http://www.cs.columbia.edu/~gusev/dash/
    to find clusters of sequences from IBD segments
    '''
    exe=os.path.join(cwd,'bin','dash_cc')
    cmds=['cat output/germline%s.match'%chrom_name,
          'cut -f 1,2,4',
          '%s output/plink%s.fam output/dash%s'%(exe,chrom_name,chrom_name)]
    try:
        output=subprocess.check_output('|'.join(cmds), shell=True)
    except subprocess.CalledProcessError as e:
        codes={1:'DASH: complete'}
        print(codes.get(e.returncode,'\nERROR CODE=%d'%e.returncode))
    except OSError:
        print('DASH executable not found at %s'%exe)

def ibd_analysis():

    # Concatenate by-chromosome GERMLINE output
    allFiles = glob.glob('output/germline**.match')
    df = pd.DataFrame()
    list = []
    for file in allFiles:
        tmp = pd.read_csv(file,delim_whitespace=True,index_col=None,header=None)
        list.append(tmp)
    df = pd.concat(list)
    df.columns=['famId1','indId1','famId2','indId2',
                'chrom','start','end','startSNP','endSNP',
                'bits','dist','unit','mismatches','homo1','homo2']

    # IBD segment lengths
    def plot_IBD_lengths(df):
        f=plt.figure('SegmentLengthIBD')
        df.dist.plot(kind='hist',bins=30,color='navy',alpha=0.2,fig=f)
        plt.xlabel('IBD segment length (cM)')

    # IBD segment locations + frequency
    def plot_IBD_map(df):
        counts=df.groupby(['chrom','start','end'])['dist'].count()
        nchrom=len(counts.index.levels[0])
        ncol = 1+nchrom/4
        nrow = int(np.ceil(float(nchrom)/ncol))
        f,axs=plt.subplots(nrow,ncol,num='SegmentMapIBD',figsize=(15,10))
        for ichrom,(chrom,ibd_counts) in enumerate(counts.groupby(level=0)):
            ax=axs[ichrom//ncol,ichrom%ncol] if nrow*ncol>1 else axs
            for idx,(rng,cnt) in enumerate(ibd_counts[chrom].iteritems()):
                ax.plot([x/1e6 for x in rng],[idx]*2,linewidth=0.1+0.6*np.log10(cnt),c='navy',alpha=0.5)
                ax.set_title('Chromosome %s'%chrom,y=0.85,x=0.25,fontsize=10)
                ax.get_yaxis().set_visible(False)
        f.set_tight_layout(True)

    # Total pairwise shared IBD lengths
    shared=df.groupby(['indId1','indId2'])['dist'].sum()
    shared.sort(ascending=False)
    def sorted_shared_idx(q):
        l=len(shared)
        return min(l-1,int(l*(1-q)))

    # IBD shared fractions
    def plot_IBD_fractions(shared):
        f=plt.figure('SharedLengthIBD')
        shared.plot(kind='hist',bins=30,color='navy',alpha=0.2,fig=f)
        plt.xlabel('Total IBD length (cM)')

    # pairwise chromosome painter
    def plot_shared_regions(df,q,ax=None,fs=12):
        idx=sorted_shared_idx(q)
        (g1,g2),v=shared.index[idx],shared[idx]
        print('q%02d genome similarity: %s and %s (IBD=%0.1f cM)'%(int(q*100),g1,g2,v))
        pair=df.groupby(['indId1','indId2']).get_group((g1,g2))
        if not ax:
            plt.figure('ChromosomePainterIBD_%s_%s'%(g1,g2))
            ax=plt.subplot(111)
        ax.set(ylim=[0,15],xlim=[-0.1,3.5],
               yticks=range(1,15),
               ylabel='chromosome',xlabel='position (MB)')
        ax.text(1.7,0.9,'\n(%s, %s)\nIBD=%0.1fcM'%(g1,g2,v),fontsize=fs)

        h=0.4
        for c,l in zip(chrom_names,chrom_lengths_Mbp):
            ax.add_patch( fpatch( (0,c-h/2), l, h,
                                   boxstyle="square,pad=0", # "round,pad=0.03"
                                   fc=(1,1,1,0),ec=(0,0,0,0.5) ) )
        for (c,s,e) in pair[['chrom','start','end']].values:
            ax.add_patch( fpatch( (s/1e6,c-h/2), (e-s)/1e6, h,
                                   boxstyle="square,pad=0",
                                   fc=(0,0,0.5,0.2),ec=(0,0,0,0.0) ) )

    def sample_shared_pairs(df,quantiles):
        L=len(quantiles)
        f,axs=plt.subplots(1,L,num='ChromosomePainterIBD',figsize=(min(18,4.5*L),4),sharex=True,sharey=True)
        for i,q in enumerate(quantiles):
            ax=axs[i] if L>1 else axs
            plot_shared_regions(df,q,ax,fs=10)
        f.set_tight_layout(True)

    # IBD fraction network
    def plot_relation_network(shared):
        f=plt.figure('RelationNetworkIBD',facecolor='w',figsize=(6,6))
        ax=f.add_axes([0,0,1,1],aspect='equal',frameon=False,xticks=[],yticks=[])
        ax.set_title('IBD network',y=0.95,color=(0.1,0.1,0.1))
        G=nx.Graph()
        ibd_spring_scale=5e-4
        for pair,ibd in shared.iteritems():
            G.add_edge(*pair,ibd=ibd,weight=ibd*ibd_spring_scale)
        gene_pop_file='output/genome_population_map.csv'
        nc=[]
        if os.path.exists(gene_pop_file):
            gene_pop_map=pd.read_csv(gene_pop_file,header=None,index_col=0,names=['gid','pid'])
            for n in G.nodes_iter():
                genome_id=int(n.replace('g',''))
                if genome_id in gene_pop_map.index:
                    nc.append(gene_pop_map.loc[genome_id]['pid'])
                else:
                    print('WARNING: no transmissions of genome id=%d'%genome_id)
                    nc.append(0)
        pos=nx.spring_layout(G,weight='weight')
        nx.draw_networkx_edges(G,pos,edge_width=1e-3,edge_color='gray',alpha=0.05)
        nc=nc if nc else 'black'
        nx.draw_networkx_nodes(G,pos,node_color=nc,node_size=60,alpha=0.4,cmap=cmap)

    # Call plotting functions
    plot_IBD_lengths(df)
    plot_IBD_map(df)
    plot_IBD_fractions(shared)
    plot_shared_regions(df,0.8)
    sample_shared_pairs(df,quantiles=[0.03,0.2,0.5,0.8,0.97])
    plot_relation_network(shared)

def cluster_analysis():
    pass

# import scipy.io
# def MATLAB_export(file='simulations/GenomeReport.npz'):
#     d=load_npz(file)
#     scipy.io.savemat('GenomeReport.mat',d)

def load_npz(file):
    try:
        with np.load(file) as data:
            A = data['genomes']
            header=data['header']
    except IOError as e:
        sys.exit(e)
    return {'genomes':A,'header':header}

def genome_analysis(file='simulations/GenomeReport.npz',reformat=True,sample=100):
    '''
    Analysis of the GenomeReport output
    '''
    d=load_npz(file)
    genomes=pd.DataFrame(d['genomes'],columns=d['header'])

    for chrom_name in chrom_names:
        if reformat:
            print('Chromosome %s:'%chrom_name)
            chrom=genomes.filter(regex='\w\.%s\.\w'%chrom_name)
            plink_format(chrom.iloc[::sample,:],chrom_name) # default: 100x downsampling
            ibd_finder(chrom_name)
            #cluster_finder(chrom_name)

    ibd_analysis()
    #cluster_analysis()

#if __name__ == '__main__':
    #MATLAB_export('../../examples/simulations/GenomeReport.npz')
    #genome_analysis('../../examples/simulations/GenomeReport.npz',
    #                reformat=True)
    #plt.show()
