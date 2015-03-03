import json
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt

with open('simulations/report.json') as reportfile:
    report=json.loads(reportfile.read())

n_infections=defaultdict(list)

sort_by_day=lambda x:int(x[0])
for tstep,data in sorted(report['data'].items(),key=sort_by_day):
    for pid,infections in data.items():
        n_infections[pid].append(len(infections))

df_infections=pd.DataFrame(n_infections)
#print(df_infections)

df_infections.plot()
plt.ylim([0,df_infections.values.max()])
plt.show()
