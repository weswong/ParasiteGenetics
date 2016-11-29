import os
import csv
import logging
log = logging.getLogger(__name__)

import numpy as np
import h5py

from .. import genome as gn
from ..infection import Transmission

class Listener:
    def notify(self,*args): pass
    def write(self, working_directory): pass

class TransmissionGeneticsReport(Listener):
    '''
    A helper class to record detailed report on genetics
    of infection transmission events to output file
    '''

    def __init__(self, parent, report_filename='TransmissionGeneticsReport.csv'):
        self.report_filename=report_filename
        self.parent=parent
        self.event='infection.transmit'
        self.data=[]

    def notify(self,*args):
        try:
            transmission=args[0]
            assert isinstance(transmission,list)
            for t in transmission:
                assert isinstance(t,Transmission)
                self.data.append(t.to_tuple())
        except:
            raise Exception('Expected list of Transmission objects as first argument.')                

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['day','pid','iid','iidParent','gidParent1','gidParent2','gid','type'])
            for r in self.data:
                writer.writerow(r)

class GenomeReport(Listener):
    '''
    A helper class to record map of Genome.id to Genome.genome
    '''

    def __init__(self, parent, report_filename='GenomeReport'):
        self.report_filename=report_filename
        self.parent=parent
        self.event='genome.init'
        self.header=gn.Genome.SNP_names
        self.data=[]

    def notify(self,*args):
        try:
            g=args[0]
            assert isinstance(g,gn.Genome)
        except:
            raise Exception('Expected Genome object as first argument.')
        barcode=g.barcode()
        self.data.append(barcode)

    # TODO: performance may be improved by doing one or both of:
    #       - store genomes instead of barcodes on notify,
    #         slicing on SNP positions in one matrix operation
    #       - pre-allocate matrix to store genomes/barcodes;
    #         trim to size before writing, but save list-to-array time

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        A=np.array(self.data)        
        with h5py.File(filename + '.h5', 'w') as hf:
            hf.create_dataset('header', data = self.header)
            hf.create_dataset('genomes', data = A)
