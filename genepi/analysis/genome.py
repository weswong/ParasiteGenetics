import numpy as np
import pandas as pd

def genome_analysis(file='simulations/GenomeReport.npz'):
    '''
    Analysis of the GenomeReport output
    '''
    with np.load(file) as data:
        A = data['genomes']
        header=data['header']
    genomes=pd.DataFrame(A,columns=header)
    print(genomes.head())

if __name__ == '__main__':
    genome_analysis('../../examples/simulations/GenomeReport.npz')
