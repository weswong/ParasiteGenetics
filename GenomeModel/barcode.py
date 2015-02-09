import genome as gn

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

allele_freqs=[ 0.7567, 0.5334, 0.9282, 0.4096,
               0.2809, 0.1045, 0.4659, 0.6340,
               0.7400, 0.1839, 0.0155, 0.1021,
               0.9507, 0.9131, 0.3911, 0.8800,
               0.8641, 0.8476, 0.1940, 0.5053,
               0.6744, 0.8951, 0.5995, 0.9217 ]

def positions_from_txt_table(filename,allele_freqs=[]):
    '''
    Read SNP positions from file in following format:
    CHR    POS
    Pf3D7_01_v3    130339
    '''
    log.info('Reading SNPs from file: %s', filename)
    SNPs=[]
    with open(filename) as f:
        for content in f.readlines()[1:]:
            CHR,POS = content.split()
            snp=gn.SNP(chrom=int(CHR.split('_')[1]),
                       pos=int(POS))
            SNPs.append(snp)

    if allele_freqs:
        if len(allele_freqs) != len(SNPs):
            raise Exception('Incompatible lengths of \
                             SNP positions and allele frequencies')
        for i,s in enumerate(SNPs):
            s.allele_freq=allele_freqs[i]

    return SNPs

SNPs=positions_from_txt_table('barcode_loci.txt',allele_freqs)
