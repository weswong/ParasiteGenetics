import json
from collections import Counter,defaultdict
import itertools
import datetime

import pandas as pd
import matplotlib.pyplot as plt

with open('simulations/report.json') as reportfile:
    report=json.loads(reportfile.read())

tsteps=[]
firstday=datetime.date(2000,1,1)

n_infections=defaultdict(list)
n_genomes=defaultdict(list)
COI={pid:{1:[],2:[],3:[],'4+':[]} for pid in report['populations']}

def COI_counts(counter):
    def COI_category(coi):
        return coi if coi<4 else '4+'
    category_counts={1:0,2:0,3:0,'4+':0}
    for k,v in counter.items():
        category_counts[COI_category(k)]+=v
    return category_counts

sort_by_day=lambda x:int(x[0])
for tstep,data in sorted(report['genomes'].items(),key=sort_by_day):
    tsteps.append(int(tstep))
    for pid,infections in data.items():
        coi_counter=Counter([len(set(i)) for i in infections])
        category_counts=COI_counts(coi_counter)
        for cat,count in category_counts.items():
            COI[pid][cat].append(count)
        n_infections[pid].append(len(infections))
        all_genomes=[i for i in itertools.chain.from_iterable(infections)]
        n_genomes[pid].append(len(set(all_genomes)))

dates=[firstday+datetime.timedelta(days=t) for t in tsteps]

f,axs = plt.subplots(4,1,figsize=(8,8),sharex=True)

df_humans=pd.DataFrame(report['n_humans'],index=dates)
#print(df_humans.head())
ax=df_humans.plot(ax=axs[0])
ax.set_ylim([df_humans.values.min(),df_humans.values.max()])
ax.set_title('Total # of humans')

df_infections=pd.DataFrame(n_infections,index=dates)
#print(df_infections.head())
ax=df_infections.plot(ax=axs[1])
ax.set_ylim([0,df_infections.values.max()])
ax.set_title('# of infected individuals')

df_genomes=pd.DataFrame(n_genomes,index=dates)
#print(df_genomes.head())
ax=df_genomes.plot(ax=axs[2])
ax.set_ylim([0,df_genomes.values.max()])
ax.set_title('# of unique genomes')

df_COI={pid:pd.DataFrame(d,index=dates) for pid,d in COI.items()}
multiple={}
for pid,df in df_COI.items():
    multiple[pid] = df.loc[:,2:"4+"].sum(axis=1).div(df.sum(axis=1), axis=0)
df_multiple=pd.DataFrame(multiple,index=dates)
#print(df_multiple.head())
ax=df_multiple.plot(ax=axs[3])
ax.set_ylim([0,1])
ax.set_title('Fraction of multiple infections')

plt.tight_layout()
plt.show()
