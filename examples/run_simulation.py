import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

import genepi.genome as gn
import genepi.simulation as sim
import genepi.report.report as report

def init_genome():
    #gn.initialize_from('barcode')
    gn.initialize_from('sequence',bin_size=1000,min_allele_freq=0.03)

def init_demographics():
    sim.Demographics.initialize_from('two_node')

def run_simulation():
    s=sim.Simulation()
    s.add_report(report.PopulationInfectionReport)
    s.add_listener(report.TransmissionGeneticsReport)
    s.add_listener(report.GenomeReport)
    s.populate_from_demographics()
    s.run()

if __name__ == '__main__':
    init_genome()
    init_demographics()
    run_simulation()
