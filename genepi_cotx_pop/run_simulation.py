import logging
import sys
import json
import shutil
log = logging.getLogger('genepi_cotx_pop')
logging.basicConfig(format='%(message)s')
log.setLevel(logging.INFO)

import genepi_cotx_pop.genome as gn
import genepi_cotx_pop.simulation as sim
from genepi_cotx_pop.report.report_v2 import PopulationInfectionReport
from genepi_cotx_pop.report.listener import TransmissionGeneticsReport,GenomeReport

from analyze_reports import run_analyses

def init_genome():
    #gn.initialize_from('barcode')
    gn.initialize_from('sequence', bin_size=500, min_allele_freq=0.03) #binsize chosen to enforce barcode frequencies

def populate_demographics(s):
    s.populate_from_demographics(node_base_name)
    #s.populate_from_demographics('multi_node', N=10,V =(0.18,3e-3,8), M=4e-4) #1e-4, 5e-4
    #s.populate_from_demographics('grid_node', L=3, V=(0.18,3e-3,8), M=2e-4)

def add_reports(s):
    s.add_reports(PopulationInfectionReport)
    s.add_listeners(TransmissionGeneticsReport,GenomeReport)

def run_simulation():
    init_genome()
    s=sim.Simulation(sim_duration=365* 30,
                     n_oocyst_params=n_oocyst_params, 
                     n_hepatocyte_params = n_hepatocyte_params,
                     infection_duration_params = infection_duration_params,
                     working_dir=working_dir+str(iter))
    add_reports(s)
    populate_demographics(s)
    s.run()
    write_params(s, working_dir)
    
def write_params(s,working_dir, param_report='SimulationParameters.json', demog_report='SimulationDemog.json'):
    param_filename =  working_dir + param_report
    params_dict = s.params.__dict__
    params_dict['demog_file'] = node_base_name
    with open(param_filename,'w') as outfile:
        json.dump(params_dict, outfile, sort_keys=True)


if __name__ == '__main__':
    node_base_name = sys.argv[1]
    iter = sys.argv[2]
    n_oocyst_params = (2.5, 0.65)
    n_hepatocyte_params = (1.6, 0.8)
    infection_duration_params = (5.13, 0.8)
    working_dir = node_base_name + '/'
    run_simulation()
