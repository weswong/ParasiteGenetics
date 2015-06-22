import unittest
from collections import Counter

from test_utils import binom_interval
import genepi.migration as mig
import genepi.population as pop
import genepi.simulation as sim
from genepi.demog.demog import annual_cycle

class TestMigration(unittest.TestCase):

    def setUp(self):
        self.rates = {'d1':0.1,'d2':0.03}
        print(self.rates)
        self.mig_info = mig.MigrationInfo(self.rates)

    def migration_destination_test(self):
        confidence_interval = 0.99
        n_humans = 100.
        dt = 2
        destinations = self.mig_info.destinations_in_timestep(n_humans, dt)
        migrants = Counter(destinations)
        print(migrants)
        for d,m in migrants.items():
            lower, upper = binom_interval(self.rates[d]*n_humans*dt, n_humans, confidence_interval)
            print('%s: %f < %f < %f ?' % (d, lower, m/float(n_humans), upper))
            self.assertTrue(lower <= m/n_humans <= upper)

class TestTwoNodeMigration(unittest.TestCase):

    migration_rate = 0.02

    def setUp(self):
        self.simulation = sim.Simulation(sim_duration = 21*6, sim_tstep = 21)

        demog = {
          'Desert' : {
            'n_humans' : 10,
            'n_infections' : 0,
            'vectorial_capacity_fn' : annual_cycle(0,coeff=0),
            'migration_rates' : {'Bog': self.migration_rate},
          },
          'Bog' : {
            'n_humans' : 10,
            'n_infections' : 2,
            'vectorial_capacity_fn' : annual_cycle(0.1,coeff=0),
            'migration_rates' : {'Desert': self.migration_rate},
          }
        }

        self.simulation.populations = { k:pop.Population(k,self.simulation,**v) \
                                        for k,v in demog.items() }

    def assertions_per_tstep(self):
        self.simulation.update()
        self.assertEqual(20, sum([p.n_humans() for p in self.simulation.populations.values()]))

    def cohort_individual_test(self):
        params = self.simulation.params
        for t in range(params.sim_duration/params.sim_tstep):
            self.assertions_per_tstep()

class TestTwoNodeWithoutMigration(TestTwoNodeMigration):

    migration_rate = 0

    def assertions_per_tstep(self):
        for p in self.simulation.populations.values():
            self.assertEqual(10, p.n_humans())
        self.assertEqual(0, self.simulation.populations['Desert'].n_infecteds())

if __name__ == '__main__':
    unittest.main()
