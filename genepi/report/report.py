import os
import csv
import json
from collections import defaultdict

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import numpy as np

from .. import genome as gn
from .. import simulation as sim

class PopulationInfectionReport:
    '''
    A helper class to record high-level metrics of simulation state
    of infections in populations to output file
    '''

    def __init__(self, parent, report_filename='PopulationInfectionReport.json'):
        self.report_filename=report_filename
        self.parent=parent
        self.data={# TODO: demographics instance owned by Simulation?
                   'populations':sim.Demographics.populations.keys(),
                   'tsteps':[],
                   'n_humans':defaultdict(list),
                   'f_infected':defaultdict(list),
                   'f_polygenomic':defaultdict(list)}

    def update(self):
        log.debug('Report updated at t=%d',self.parent.day)
        self.data['tsteps'].append(self.parent.day)
        for pid,p in self.parent.populations.items():
            n_humans=p.n_humans()
            self.data['n_humans'][pid].append(n_humans)
            n_infecteds=p.n_infecteds()
            f_infected=float(n_infecteds)/n_humans if n_humans else 0
            self.data['f_infected'][pid].append(f_infected)
            f_poly=float(p.n_polygenomic())/n_infecteds if n_infecteds else 0
            self.data['f_polygenomic'][pid].append(f_poly)

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename,'w') as outfile:
            json.dump(self.data, outfile, sort_keys=True)

class TransmissionGeneticsReport:
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

class GenomeReport:
    '''
    A helper class to record map of Genome.id to Genome.genome
    '''

    def __init__(self, parent, report_filename='GenomeReport.csv'):
        self.report_filename=report_filename
        self.parent=parent
        self.event='genome.init'
        filename=os.path.join(sim.Params.working_dir, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        self.csvfile = open(filename, 'wb')
        self.writer = csv.writer(self.csvfile)
        #self.writer.writerow(['id']+gn.Genome.SNP_names) # header
        self.writer.writerow(gn.Genome.SNP_names) # header

    @profile
    def notify(self,*args):
        try:
            g=args[0]
        except:
            raise Exception('Expected Genome object as first argument.')
        barcode=g.barcode()
        self.writer.writerow(barcode)

    def write(self, working_directory):
        self.csvfile.close()
