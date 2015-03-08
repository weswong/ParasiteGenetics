import os
import csv

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

from genepi import genome as gn

def positions_from_csv(filename):
    '''
    Read SNP positions from a CSV file in following format:
    chrom,pos,freq
    1,1110,0.02631
    1,2271,0.03125
    1,2402,0.03125
    1,4559,0.02777
    '''
    log.info('Reading SNPs from file: %s', filename)
    SNPs=[]
    with open(os.path.join(os.path.dirname(__file__),filename)) as f:
        reader=csv.reader(f)
        header=next(reader,None)
        for chrom,pos,freq in reader:
            snp=gn.SNP(int(chrom),int(pos),float(freq))
            SNPs.append(snp)

    return SNPs

SNPs=positions_from_csv('sequence_loci.csv')
