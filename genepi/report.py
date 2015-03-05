import os
import json
from collections import defaultdict

import logging
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

import genome as gn
import simulation as sim # TODO: demographics instance owned by Simulation?

class SimulationReport():
    '''
    A helper class to keep track of simulation states to be reported
    '''

    def __init__(self, parent, report_filename='report.json'):
        self.report_filename=report_filename
        self.parent=parent
        self.data={'SNPs':[s.to_tuple() for s in gn.Genome.SNPs],
                   'populations':sim.Demographics.populations.keys(),
                   'genomes':{},
                   'n_humans':defaultdict(list)}

    def update(self):
        log.debug('Report updated at t=%d',self.parent.day)
        today=self.data['genomes'][self.parent.day]=defaultdict(list)
        for pid,iid,i in self.parent.iterate_infections():
            if i:
                today[pid].append([g.barcode_as_long() for g in i.genomes])
            else:
                today[pid].append([])
        for pid,p in self.parent.populations.items():
            self.data['n_humans'][pid].append(p.n_humans())

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename,'w') as outfile:
            json.dump(self.data, outfile, sort_keys=True)
