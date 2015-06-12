from functools import partial

import logging
log = logging.getLogger('genepi')
logging.basicConfig(format='%(message)s')
log.setLevel(logging.INFO)

import pandas as pd
import matplotlib.pyplot as plt

import genepi.genome as gn
import genepi.simulation as sim
from genepi.report.report import PopulationInfectionReport
from genepi.report.listener import TransmissionGeneticsReport,GenomeReport
from genepi.analysis.genome import genome_analysis, load_npz
from genepi.analysis.population import population_analysis
from genepi.analysis.transmission import transmission_analysis
from genepi.event.drug import resistant_sites, set_resistance, DrugTreatment

def init_genome():
    #gn.initialize_from('barcode')
    gn.initialize_from('sequence',bin_size=1000,min_allele_freq=0.03)

    #AM_resistance
    chrom,pos = resistant_sites['AM']
    gn.add_locus(chrom,pos,'Pf.%d.%d'%(chrom,pos),fitness=0.95)

def populate_demographics(s):
    #s.populate_from_demographics('single_node')
    s.populate_from_demographics('multi_node',N=10,V=(0.18,4e-3,8),M=5e-4)

def add_reports(s):
    s.add_reports(PopulationInfectionReport)

    treatment_by_population = lambda t: {'Population1': {
                                             'fraction': max(0.8,0.3+0.1*t/365),
                                             'clearance': 0.95,
                                             'resistant_clearance': 0.5 }}

    s.add_listeners(TransmissionGeneticsReport,
                    GenomeReport,
                    partial(DrugTreatment,
                            drug='AM',
                            treatment_by_population=treatment_by_population))

def run_simulation():
    init_genome()
    s=sim.Simulation(sim_duration=365*10) #365*20
    add_reports(s)
    populate_demographics(s)

    tstep=s.params.sim_tstep
    #for t in range(28*tstep,700,tstep):
    for i,t in enumerate(range(28*tstep,35*tstep,tstep)):
        resistant_genome=gn.Genome.from_allele_freq(mod_fns=[partial(set_resistance,drug='AM',allele=i+1)])
        mutation_fn = lambda s: s.populations.values()[0].add_infection_from_genomes([resistant_genome])
        s.add_event(t, mutation_fn)

    s.run()

def resistance_analysis(res_file='simulations/GenomeReport.npz',
                        tx_file='simulations/TransmissionGeneticsReport.csv'):

    d=load_npz(res_file)
    genomes=pd.DataFrame(d['genomes'],columns=d['header'])
    resistance_locus=genomes.keys()[-1]
    print('Resistance',resistance_locus)
    genome_resistance=genomes[resistance_locus]
    resistant_genomes=genome_resistance[genome_resistance>0].index.values
    print('# resistant genomes observed',len(resistant_genomes))

    tx=pd.read_csv(tx_file)
    resistant_tx = tx[tx.gid.isin(resistant_genomes)]

    resistant_incidence=resistant_tx.groupby(['pid','day'])['iid'].nunique().unstack('pid')
    resistant_incidence.plot(title='AM resistant incidence')

if __name__ == '__main__':
    run_simulation()
    population_analysis()
    transmission_analysis()
    genome_analysis()
    resistance_analysis()
    plt.show()
