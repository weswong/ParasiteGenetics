import os
import json
from collections import defaultdict

import logging
log = logging.getLogger(__name__)

class Report:
    def update(self): pass
    def write(self, working_directory): pass

class PopulationInfectionReport(Report):
    '''
    A helper class to record high-level metrics of simulation state
    of infections in populations to output file
    '''

    def __init__(self, parent, report_filename='PopulationInfectionReport.json'):
        self.report_filename=report_filename
        self.parent=parent
        self.data={'tsteps':[],
                   'n_humans':defaultdict(list),
                   'f_infected':defaultdict(list),
                   'f_polygenomic':defaultdict(list),
                   'vectorial_capacity':defaultdict(list),
                   'coi_distribution': defaultdict(list)}

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
            
            vectorial_capacity = p.vectorial_capacity()
            self.data['vectorial_capacity'][pid].append(vectorial_capacity)
            
            coi_distribution = p.coi_distribution()
            self.data['coi_distribution'][pid].append(coi_distribution)

    def write(self, working_directory):
        self.data['populations']=self.parent.populations.keys()
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename,'w') as outfile/:
            json.dump(self.data, outfile, sort_keys=True)
