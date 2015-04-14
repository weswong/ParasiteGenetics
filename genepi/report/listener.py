import os
import csv
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import numpy as np

from .. import genome as gn

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
        except:
            raise Exception('Expected Transmission object as first argument.')
        self.data.append(transmission.to_tuple())

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['day','pid','iid','iidParent','gidParent1','gidParent2','gid'])
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
        except:
            raise Exception('Expected Genome object as first argument.')
        barcode=g.barcode()
        self.data.append(barcode)

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        A=np.array(self.data)
        np.savez(filename,genomes=A,header=self.header)
