from functools import partial
import datetime

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
import genepi.event.drug as drug

def init_genome():
    gn.initialize_from('sequence', bin_size=1000, min_allele_freq=0.03)
    chrom, pos = drug.resistant_sites['AM']
    gn.add_locus(chrom, pos, 'Pf.%d.%d' % (chrom,pos), fitness=0.95, freq=0)

def populate_demographics(s):
    s.populate_from_demographics('grid_node', L=3, V=(0.18,3e-3,8), M=5e-3)

def add_reports(s):
    s.add_reports(PopulationInfectionReport)
    s.add_listeners(TransmissionGeneticsReport, GenomeReport)

def add_treatment(s):
    treatment_by_population = lambda t: {p: {'fraction': min(0.6, 0.05 + 0.05*t/365),
                                             'clearance': lambda g: 0.05 if drug.is_resistant(g, 'AM') else 0.9} for p in s.populations.keys()}
    s.add_listeners(partial(drug.DrugTreatment,
                            treatment_by_population=treatment_by_population))

def run_simulation():
    init_genome()
    s=sim.Simulation(sim_duration=365*10)
    add_reports(s)
    populate_demographics(s)
    add_treatment(s)

    for i,t in enumerate(range(0, 365, s.params.sim_tstep)):
        mf = partial(drug.set_resistance, drug='AM') # (allele=i+1) but GERMLINE wants bi-allelic instead of marking individual emergences with 1,2,3,etc.
        resistant_genome=gn.Genome.from_allele_freq(mod_fns = [mf])
        mutation_fn = lambda s: s.populations['Population #0'].add_infection_from_genomes([resistant_genome])
        s.add_event(t, mutation_fn)

    s.run()

def resistance_analysis(res_file='simulations/GenomeReport.npz',
                        tx_file='simulations/TransmissionGeneticsReport.csv'):

    d=load_npz(res_file)
    genomes=pd.DataFrame(d['genomes'],columns=d['header'])
    resistance_locus='Pf.%d.%d' % drug.resistant_sites['AM']
    print('Resistance',resistance_locus)
    genome_resistance=genomes[resistance_locus]
    resistant_genomes=genome_resistance[genome_resistance>0].index.values
    print('# resistant genomes observed',len(resistant_genomes))

    tx=pd.read_csv(tx_file)
    resistant_tx = tx[tx.gid.isin(resistant_genomes)]

    resistant_incidence=resistant_tx.groupby(['pid','day'])['iid'].nunique().unstack('pid').fillna(0)

    #print(resistant_incidence)

    firstday=datetime.date(2000,1,1)
    dates=[firstday+datetime.timedelta(days=t) for t in resistant_incidence.index]
    resistant_incidence.index = dates
    resistant_incidence.plot(title='AM resistant incidence')

if __name__ == '__main__':
    run_simulation()
    population_analysis()
    transmission_analysis()
    genome_analysis(sample=1)
    resistance_analysis()
    plt.show()
