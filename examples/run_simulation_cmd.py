import logging
log = logging.getLogger('genepi')
logging.basicConfig(format='%(message)s')
log.setLevel(logging.INFO)

import sys
import genepi.genome as gn
import genepi.simulation as sim
from genepi.report.report import PopulationInfectionReport
from genepi.report.listener import TransmissionGeneticsReport,GenomeReport

from analyze_reports import run_analyses

def init_genome():
    #gn.initialize_from('barcode')
    gn.initialize_from('sequence', bin_size=1000, min_allele_freq=0.03)

def populate_demographics(s):
    s.populate_from_demographics('single_node_{number}'.format(number=demog_file_number))
    #s.populate_from_demographics('multi_node', N=10,V =(0.18,3e-3,8), M=4e-4) #1e-4, 5e-4
    #s.populate_from_demographics('grid_node', L=3, V=(0.18,3e-3,8), M=2e-4)

def add_reports(s):
    s.add_reports(PopulationInfectionReport)
    s.add_listeners(TransmissionGeneticsReport,GenomeReport)

def run_simulation():
    init_genome()
    s=sim.Simulation(sim_duration=365*25, working_dir = 'simulations_{number}'.format(number=demog_file_number))
    add_reports(s)
    populate_demographics(s)
    s.run()

if __name__ == '__main__':
    demog_file_number = sys.argv[1]
    run_simulation()
    #run_analyses()
