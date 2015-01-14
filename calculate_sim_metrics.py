from parse_sim import BarcodeSimulationParser

class BarcodeMetricCalculator():

    def __init__(self, work_dir='simulations', idx=0):
        self.parser=BarcodeSimulationParser(work_dir, idx)
        self.samples=self.parser.sample_monogenomic_barcodes()
        (self.repeats,self.persists)=self.parser.get_repeated_and_persistent_strains(self.samples)

        self.sim={}
        self.sim_err={}
        self.calculate_mixed_fraction()
        self.calculate_persistence()
        self.calculate_repeats()

    def __str__(self):
        s=[ 'Mixed Fraction: [%s]' % ', '.join('%.3f' % val for val in self.sim['mixed_infection_fraction']),
            'Created: ' + str(self.sim['created_persists']),
            'Destroyed: ' + str(self.sim['destroyed_persists']),
            'Persist Years: ' + str(self.sim['persist_years']),
            'Unique Fraction: [%s]' % ', '.join('%.3f' % val for val in self.sim['unique_fraction']),
            'Size-2 Repeats: ' + str(self.sim['repeats_equal_two']),
            'Greater-than-2 Repeat per Repeat: [%s]' % ', '.join('%.3f' % val for val in self.sim['greater_two_per_repeat']) ]
        return '\n'.join(s)

    def calculate_mixed_fraction(self):
        mf=self.parser.fraction_mixed[0][::365]      # one ParasitePopulation, one value per year
        self.sim['mixed_infection_fraction']=mf[-8:] # last eight years map to 2006-2013 in Thies

        # Try to account for variation within malaria season and potential bias
        # of clinical cases towards first new infection of season.
        mixed_fraction_in_season=[]
        for yr in range(self.parser.n_years)[-8:]:
            season=self.parser.fraction_mixed[0][yr*365-180:yr*365]
            mixed_fraction_in_season.append((min(season),max(season)))

        self.sim['mixed_infection_fraction'] = [0.5*(mn+mx) for mn,mx in mixed_fraction_in_season]
        self.sim_err['mixed_infection_fraction'] = [0.5*(mx-mn) for mn,mx in mixed_fraction_in_season]

    def calculate_persistence(self):
        persist_counts = [0,0,0,0] # 2,3,4,5+ year spans
        created,destroyed=0,0
        for persists in self.persists.values():
            if persists[0]==0:
                created += 1
            if persists[-1]==0:
                destroyed += 1
            nonzero_indices=[i for i, e in enumerate(persists) if e > 0]
            persist_range = nonzero_indices[-1]-nonzero_indices[0]+1
            persist_counts[min(persist_range,5)-2] += 1 # count all >= 5 years together

        self.sim['persist_years']=persist_counts
        self.sim['created_persists']=created
        self.sim['destroyed_persists']=destroyed

    def calculate_repeats(self):
        self.sim['unique_fraction']=[]
        self.sim['repeats_equal_two']=[]
        self.sim['greater_two_per_repeat']=[]
        for repeats_in_year in self.repeats: # loop over years
            repeat_counts = repeats_in_year.values()
            n_samples = sum(repeat_counts)
            n_unique = sum([r==1 for r in repeat_counts])
            n_two = sum([r==2 for r in repeat_counts])
            n_greater = sum([r>2 for r in repeat_counts])
            self.sim['unique_fraction'].append(n_unique/float(n_samples) if n_samples>0 else 0)
            self.sim['repeats_equal_two'].append(n_two)
            self.sim['greater_two_per_repeat'].append(n_greater/float(n_two+n_greater) if (n_two+n_greater)>0 else 0)

if __name__ == '__main__':
    c=BarcodeMetricCalculator()
    print(c)
