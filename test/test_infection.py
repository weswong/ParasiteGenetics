import unittest

import numpy as np
import matplotlib.pyplot as plt
from nose.tools import nottest

from test_utils import binom_interval
import genepi.genome as gn
import genepi.infection as inf
import genepi.utils as utils

@nottest
def plot_chokepoints_test():
    import pandas as pd
    import seaborn as sns

    o = [inf.sample_n_oocysts() for _ in range(2000)]
    h = [inf.sample_n_hepatocytes() for _ in range(2000)]
    data=pd.DataFrame({'oocysts':o, 'hepatocytes':h})

    sns.set(style='ticks')

    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.hist(o, bins=np.arange(-0.5, 49.5), alpha=0.3)
    ax1.set_title('Oocysts')
    ax2.hist(h, bins=np.arange(-0.5, 49.5), alpha=0.3)
    ax2.set_title('Hepatocytes')
    ax2.set_xlim([-1, 40])

    sns.jointplot(x='oocysts', y='hepatocytes', data=data,
                  kind='scatter', stat_func=None, dropna=False, color='navy',
                  marginal_kws={'bins':np.arange(-0.5,49.5)},
                  joint_kws={'alpha':0.1},
                  xlim=[-1,40], ylim=[-1,40])


class TestInfectionTransmit(unittest.TestCase):

    nRepeats = 5
    nSigma = 3

    def setUp(self):
        gn.initialize_from('barcode', bin_size = 1e7)

    def test_transmit(self):
        g1 = gn.Genome.from_reference()
        g2 = gn.Genome.from_barcode([1]*gn.num_SNPs())

        i1=inf.Infection(None, [g1])
        tt = i1.transmit()
        self.assertEqual(len(tt), 1)
        self.assertEqual(tt[0].genome.id, g1.id)

        i2=inf.Infection(None, [g1, g2])
        for _ in range(self.nRepeats):
            tt = i2.transmit()
            for tx in tt:
                print(tx.parentGenomeIds)
                print(tx.parentInfection)
                self.assertEqual(tx.parentInfection.id, i2.id)
                self.assertTrue(all([g in [g1.id, g2.id] for g in tx.parentGenomeIds]))

    def test_oocyst_sample(self, n_hep=5, n_ooc=3):
        for _ in range(self.nRepeats):
            products = inf.sample_oocyst_products(n_hep,n_ooc)
            self.assertTrue(1 <= len(products) <= n_ooc)
            self.assertTrue(1 <= sum(products) <= n_hep)

    def generator_test(self):
        G=inf.infectious_generator()
        G.send(None)
        dt=21 # less than one incubation period
        for i in range(self.nRepeats):
            G.send(dt)
            n=next(G)
            print(n)
            if i==0:
                self.assertEqual(n, 0) # still incubating after first dt
            else:
                self.assertGreater(n, 0)

    def infection_duration_test(self):
        d = []
        for _ in range(self.nRepeats):
            d.append(np.log(inf.get_infection_duration()))
        self.assertAlmostEqual(np.mean(d), 5.13, delta = self.nSigma * 5.13 / np.sqrt(self.nRepeats))
        self.assertAlmostEqual(np.std(d), 0.8, delta = self.nSigma * 0.8 / np.sqrt(self.nRepeats))

    def reinfection_test(self):
        g1 = gn.Genome.from_reference()
        g2 = gn.Genome.from_barcode([1]*gn.num_SNPs())
        i = inf.Infection(None, [g1, g2])
        t, dt = 0, 21
        self.assertEqual(next(i.infectiousness), 0)
        for _ in range(self.nRepeats):
            txs = i.update(dt, 1)
            for tx in txs:
                nstrains = i.n_strains()
                genomes = [g.genome for g in tx]
                i.merge_infection(genomes)
                self.assertGreater(next(i.infectiousness), 0) # reinitialized post-incubation
                self.assertGreaterEqual(nstrains + len(genomes), i.n_strains())
                self.assertEqual(i.n_strains(), len(set([g.display_barcode() for g in i.genomes])))
                break # just one transmission per iteration in this test suite
            print('t=%d\n%s' % (t, str(i)))
            t += dt
        self.assertGreater(i.n_strains(), 2)

    def test_gametocyte_sampling(self):
        fitcost = 0.2

        gn.add_locus(7, 100, 'TEST', fitness = fitcost)
        g1 = gn.Genome.from_reference()
        g2 = gn.Genome.from_barcode([1]*gn.num_SNPs())
        i = inf.Infection(None, [g1, g2])
        w = i.gametocyte_strain_cdf()

        self.assertListEqual(w, [1/(1+fitcost), 1.0])

        nSamples = 30
        confInt = 0.99
        selected = [i.select_strain(w).id == g1.id for _ in range(nSamples)]
        lower, upper = binom_interval(nSamples/(1+fitcost), nSamples, confInt)
        print(lower, upper, sum(selected))
        self.assertTrue(lower <= sum(selected)/float(nSamples) <= upper)

if __name__ == '__main__':
    plot_chokepoints_test()
    plt.show()
    unittest.main()
