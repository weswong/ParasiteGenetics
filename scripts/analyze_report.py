import json
from collections import Counter,defaultdict
import itertools
import datetime

import pandas as pd
import matplotlib.pyplot as plt

# TransmissionGeneticsReport
df=pd.read_csv('simulations/TransmissionGeneticsReport.csv',index_col=['iidParent','iid'])
print(df[['gidParent1','gidParent2','gid']][-30:])

# PopulationInfectionReport
with open('simulations/PopulationInfectionReport.json') as reportfile:
    report=json.loads(reportfile.read())

firstday=datetime.date(2000,1,1)

dates=[firstday+datetime.timedelta(days=t) for t in report['tsteps']]

f,axs = plt.subplots(3,1,figsize=(8,8),sharex=True)

channels=['n_humans','f_infected','f_polygenomic']
titles=['Total # of humans','Fraction infected','Fraction polygenomic infections']
for i,c in enumerate(channels):
    df=pd.DataFrame(report[c],index=dates)
    ax=df.plot(ax=axs[i])
    ax.set_ylim([df.values.min(),df.values.max()])
    ax.set_title(titles[i])

f.set_tight_layout(True)
plt.show()
