import os
import json
import math
import random
import numpy as np # random.negative_binomial
from collections import defaultdict
from ParasitePopulation import ParasitePopulation

def mixed_infection_prob(mean_COI):
     # 1 - P(0) - P(1) = P(>=2)
     return 1 - math.exp(-mean_COI) - mean_COI*math.exp(-mean_COI)

class SimulationParameters():
    ''' A helper class to store configurable parameters of the simulation '''

    def __init__(self, param_dict):
        for k,v in param_dict.items():
            setattr(self,k,v)

    def __str__(self):
        return 'SimulationParameters: %s' % self.__dict__

class SimulationReport():
    ''' A helper class to keep track of simulation states to be reported '''

    def __init__(self, n_populations, report_filename):
        # Annual sampling of barcodes
        self.barcode_census=defaultdict(list)
        self.report_filename=report_filename

        # Time-series of various quantities
        self.n_parasites = [[] for k in range(n_populations)]
        self.n_mixed     = [[] for k in range(n_populations)]
        self.n_barcodes  = [[] for k in range(n_populations)]
        self.p_outcross  = [[] for k in range(n_populations)]

    def update(self, populations, t):
        for i,p in enumerate(populations):
            self.n_parasites[i].append(len(p.parasites))
            self.n_mixed[i].append(sum([k.mixed_mask !=0 for k in p.parasites]))
            self.n_barcodes[i].append(len(p.barcode_counts()))
            self.p_outcross[i].append(mixed_infection_prob(p.mean_COI()))
            if t % 365 == 0:
                self.barcode_census[t//365].extend(p.infections())

    def write(self, working_directory):
        filename=os.path.join(working_directory, self.report_filename)
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        with open(filename,'w') as outfile:
            json.dump({'barcode_census':self.barcode_census,
                       'n_parasites': self.n_parasites,
                       'n_barcodes': self.n_barcodes,
                       'n_mixed': self.n_mixed,
                       'p_outcross': self.p_outcross}, outfile)

def annual_cycle(seasonal_min, seasonal_coeff, days_per_cycle=365):
    def f(t):
        if seasonal_coeff == 0:
            return 1.0
        else:
            return seasonal_min + (1.0-seasonal_min)*pow(math.cos(3.1416*t/days_per_cycle), seasonal_coeff)
    return f

class ParasiteSimulation():
    ''' A simulation of parasite population dynamics '''

    def __init__(self, params):

        self.params = SimulationParameters(params)
        self.params.n_populations = int(self.params.n_populations)

        if 'random_seed' in params:
            seed=params['random_seed']
            random.seed(seed)
            report_name = 'barcode_report_%d.json' % seed
        else:
            report_name = 'barcode_report.json'

        self.populations = [ ParasitePopulation(n_humans=self.params.n_humans//self.params.n_populations, 
                                                n_parasites=self.params.initial_parasites//self.params.n_populations, 
                                                n_unique=self.params.initial_unique//self.params.n_populations,
                                                propagation_func=self.params.propagation_function
                                                ) for k in range(self.params.n_populations) ]

        self.reports = [ SimulationReport(self.params.n_populations, 
                                          report_filename=report_name) ]

        self.working_directory = 'simulations'

    @classmethod
    def from_config_file(cls,file,seed=None):

        with open(file,'r') as infile:

            params=json.loads(infile.read())

            s = params.pop('seasonality')
            params['seasonal_forcing'] = annual_cycle(s['minimum'], s['coefficient'])

            prop_func = params.pop('propagation_function', ('CONSTANT',1))
            if prop_func[0] == 'CONSTANT':
                f=lambda:prop_func[1]
            elif prop_func[0] == 'NBINOMIAL':
                f=lambda: np.random.negative_binomial(prop_func[1],prop_func[2])
            elif prop_func[0] == 'BIMODAL':
                f=lambda: prop_func[1] if random.random() < 1.0/prop_func[1] else 0
            else:
                raise Exception('Only supporting CONSTANT, NBINOMIAL, BIMODAL propagation functions.')
            params['propagation_function'] = f

            if seed is not None:
                params['random_seed'] = seed

        return cls(params)

    @classmethod
    def from_default(cls):
        params={ 'simulation_duration': 365*15,
                 'R0_initial': 3,
                 'propagation_function': lambda: 20 if random.random() < 1.0/20 else 0,
                 'rates': { 'expire': 1/180.,
                            'import': 1/100. },
                 'seasonal_forcing': annual_cycle(seasonal_min=0.1, seasonal_coeff=6),
                 'R0_reduction_year': 5,
                 'R0_final': 1.45,
                 'R0_plateau_year': 11,
                 'n_populations': 1,
                 'initial_parasites': 1500, # total for all n_populations
                 'initial_unique': 400,
                 'n_humans': 2000 }
        return cls(params)

    def run(self):
        for t in range(self.params.simulation_duration):
            self.update(t)
        self.save()

    def update(self,t):

        if t % 365 == 0:
            print('t=%d' % t)

        for p in self.populations:
            rates = self.params.rates
            outcrossing_prob = mixed_infection_prob(p.mean_COI())
            reproduction_rate = self.params.R0_initial * rates['expire'] * self.params.seasonal_forcing(t) * max(0, 1.0 - p.mean_COI()/5.0)
            R0_reduction_scale = self.params.R0_final / self.params.R0_initial
            if t > self.params.R0_reduction_year*365 and t < self.params.R0_plateau_year*365:
                reproduction_rate *= 1 - (t-self.params.R0_reduction_year*365) * (1-R0_reduction_scale) / (self.params.R0_plateau_year*365-self.params.R0_reduction_year*365)
            elif t > self.params.R0_plateau_year*365 and t < self.params.R0_rebound_year*365:
                reproduction_rate *= R0_reduction_scale
            elif t > self.params.R0_rebound_year*365:
                reproduction_rate *= self.params.R0_rebound / self.params.R0_initial
            rates.update( { 'outcross':  outcrossing_prob * reproduction_rate,
                            'clone': (1-outcrossing_prob) * reproduction_rate } )
            p.update(rates)

        for r in self.reports:
            r.update(self.populations, t)

    def save(self):
        for r in self.reports:
            r.write(self.working_directory)
