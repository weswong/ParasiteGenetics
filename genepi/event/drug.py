import random
import logging
log = logging.getLogger(__name__)

from .. import genome as gn
from genepi.report.listener import Listener
from genepi.infection import Transmission

resistant_sites = {'AM': (13,1886271)}

def is_resistant(g,drug):
    resistant_idx = gn.snp_bin_from_chrom_pos(*resistant_sites[drug])
    return g[gn.Genome.SNP_bins[resistant_idx]]

def set_resistance(g,drug,allele=1):
    resistant_idx = gn.snp_bin_from_chrom_pos(*resistant_sites[drug])
    log.debug('Setting %s resistance at index %d',drug,resistant_idx)
    g[resistant_idx] = allele

class DrugTreatment(Listener):
    '''
    Treatment with a drug on a new-infection transmission event

    e.g. treatment_by_population = lambda t: {'Population1': {
                                                  'fraction': max(0.8,0.3+0.1*t/365),
                                                  'clearance': 0.95,
                                                  'resistant_clearance': 0.5 }}
    '''

    def __init__(self, parent, drug, treatment_by_population):
        self.parent=parent
        self.drug=drug
        self.treatment_by_population=treatment_by_population
        self.event='infection.transmit'

    def notify(self,*args):
        try:
            transmission=args[0]
            assert isinstance(transmission,list)
        except:
            raise Exception('Expected list of Transmission objects as first argument.')

        popId = transmission[0].populationId
        treatment = self.treatment_by_population(self.parent.day).get(popId,None)

        if not treatment or random.random() > treatment['fraction']:
            return

        infection = transmission[0].infection
        infection.genomes = [g for g in infection.genomes if random.random()>treatment['clearance'](g.genome)]
        if not infection.genomes:
            log.debug('New infection cleared with prompt treatment.')
            infection.infection_timer = 0
