import json
from collections import Counter,defaultdict
import itertools
import datetime
import sys

import pandas as pd
import matplotlib.pyplot as plt

try:
    import seaborn as sns
    sns.set_palette('husl',10)
    sns.set_style("white",{'grid.color': '.93'})
    sns.set_context("notebook")
except:
    print('Install seaborn package for more pleasing aesthetics.')
    pass

def population_analysis(file='simulations/PopulationInfectionReport.json'):
    '''
    Analysis of the PopulationInfectionReport output
    '''
    try:
        with open(file) as reportfile:
            report=json.loads(reportfile.read())
    except IOError as e:
        sys.exit(e)

    firstday=datetime.date(2000,1,1)

    dates=[firstday+datetime.timedelta(days=t) for t in report['tsteps']]
    print dates
    #f,axs = plt.subplots(3,1,figsize=(12,8),sharex=True)
    f,axs = plt.subplots(2,1,figsize=(12,8),sharex=True)
    
    channels=['n_humans','f_infected','f_polygenomic']
    #titles=['Total # of humans','Fraction infected', 'Fraction polygenomic infections']
    titles=['Total # of humans','Fraction infected & Fraction polygenomic infections']
    for i,c in enumerate(channels):
        df=pd.DataFrame(report[c],index=dates)
        ax=df.plot(ax=axs[i],legend=False,linewidth=1)
        ax.set_ylim([df.values.min(),df.values.max()])
        ax.set_title(titles[i])

    f.set_tight_layout(True)

if __name__ == '__main__':
    population_analysis('../../examples/simulations/PopulationInfectionReport.json')
    plt.show()
