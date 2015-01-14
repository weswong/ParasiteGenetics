import os
import json
import random
import warnings
from Parasite import Parasite
from collections import Counter

class BarcodeSimulationParser():

    def __init__(self, work_dir='simulations', idx=0):

        print('Initialize simulation parser from output directory: %s (idx=%d)' % (work_dir,idx))
        self.work_dir = work_dir
        self.sim_idx = idx

        with open(os.path.join(work_dir,'barcode_report_%d.json' % idx),'r') as infile:
            bcj = json.loads(infile.read())

        self.n_parasites = bcj['n_parasites']
        self.n_barcodes = bcj['n_barcodes']
        self.n_mixed = bcj['n_mixed']
        self.p_outcross = bcj['p_outcross']
        self.fraction_mixed = [[m/float(p) if p>0 else 0 for (m,p) in zip(my,py)] for (my,py) in zip(self.n_mixed,self.n_parasites)]
        self.census = bcj['barcode_census']
        self.n_years = len(self.census)

        self.annual_samples = [100]*self.n_years
        self.annual_samples[-8:] = [90, 138, 96, 81, 100, 112, 158, 236] # using Thies sampled numbers in last 8 years

    def sample_monogenomic_barcodes(self, years=[], seed=1234):
        if not years:
            years=range(self.n_years-8, self.n_years)
        monogenomic_samples_by_year=[]
        random.seed(seed)
        for year in years:
            annual_infections = self.census[str(year)]
            all_monogenomic = [s for s in annual_infections if not s[1]] # sampling that chooses samples from set with mixing_mask=0
            if self.annual_samples[year] > len(all_monogenomic):
                warnings.warn('Not enough infections in simulation to match data samples.')
                sampled_monogenomic = all_monogenomic
            else:
                sampled_monogenomic = random.sample(all_monogenomic, self.annual_samples[year])
            monogenomic_samples_by_year.append(sampled_monogenomic)
        return monogenomic_samples_by_year

    def get_repeated_and_persistent_strains(self, samples_by_year):
        nyears=len(samples_by_year)
        repeats_by_year = [ Counter( [ p[0] for p in samples ] ) for samples in samples_by_year ]
        strains_by_year = [ r.keys() for r in repeats_by_year ]
        all_strains = set([s for y in strains_by_year for s in y])
        strains_across_years = {s:[0]*nyears for s in all_strains}
        for i,r in enumerate(repeats_by_year):
            for k,v in r.items():
                strains_across_years[k][i] += v
        persistent_strains={k:v for k,v in strains_across_years.items() if sum([n>0 for n in v]) > 1}
        return repeats_by_year,persistent_strains

if __name__ == '__main__':
    b=BarcodeSimulationParser()
    s=b.sample_monogenomic_barcodes()
    r,p=b.get_repeated_and_persistent_strains(s)
    print(p)
