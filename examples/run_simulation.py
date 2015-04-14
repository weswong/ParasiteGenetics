import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

import genepi.genome as gn
import genepi.simulation as sim
import genepi.report.report as report

def init_genome():
    #gn.initialize_from('barcode')
    gn.initialize_from('sequence',bin_size=1000,min_allele_freq=0.03)

def populate_demographics(s):
    #s.populate_from_demographics('single_node')
    s.populate_from_demographics('multi_node',N=10,V=(0.18,1e-3,8),M=5e-4)

def add_reports(s):
    s.add_report(report.PopulationInfectionReport)
    s.add_listener(report.TransmissionGeneticsReport)
    s.add_listener(report.GenomeReport)

def run_simulation():
    init_genome()
    s=sim.Simulation(sim_duration=365*20)
    add_reports(s)
    populate_demographics(s)
    s.run()

if __name__ == '__main__':
    run_simulation()
