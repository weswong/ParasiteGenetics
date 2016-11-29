from genepi_cotx_pop.analysis.genome import genome_analysis
from genepi_cotx_pop.analysis.population import population_analysis
from genepi_cotx_pop.analysis.transmission import transmission_analysis

import matplotlib.pyplot as plt

def run_analyses():
    population_analysis()
    transmission_analysis()
    #genome_analysis()    
    plt.show()

if __name__ == '__main__':
    run_analyses()
