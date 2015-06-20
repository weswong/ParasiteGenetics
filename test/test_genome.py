import math
import unittest

import numpy as np

import genepi.genome as gn

class TestGenome(unittest.TestCase):

    nRandom = 30
    nSigma = 3

    def setUp(self):
        gn.initialize_from('barcode')

    def test_crossover_points(self):
        chrom1size = gn.Genome.chrom_breaks[1]

        firstxpoints=[]
        for _ in range(self.nRandom):
            xpoints = gn.get_crossover_points(chrom1size)
            if xpoints:
                firstxpoints.append(xpoints[0])
        lambd = gn.bp_per_morgan / gn.Genome.bin_size_bp
        k = chrom1size / lambd
        expect_trunc_expo = lambd * ((1 - (k + 1) * math.exp(-k)) / (1-math.exp(-k)))
        self.assertAlmostEqual(np.mean(firstxpoints),
                               expect_trunc_expo,
                               delta = self.nSigma*expect_trunc_expo/math.sqrt(self.nRandom))

    def test_meiosis(self):
        in1 = gn.Genome.from_reference()
        in2 = gn.Genome.from_barcode([1]*gn.num_SNPs())
        print(in1.display_barcode())
        print(in2.display_barcode())
        genomes = gn.meiosis(in1, in2)
        print('\n'.join([g.display_barcode() for g in genomes]))
        self.assertListEqual(np.sum([g.barcode() for g in genomes], axis=0).tolist(), [2]*gn.num_SNPs())
        self.assertRaises(IndexError, gn.meiosis, in1, in2, N=5)

    def test_gametocyte_products(self):
        g1 = gn.Genome.from_reference()
        g2 = gn.Genome.from_barcode([1]*gn.num_SNPs())

        tt = gn.distinct_sporozoites_from([(g1,g1)], n_products=[3])
        self.assertEqual(len(tt), 1)
        self.assertEqual(tt[0].genome.id, g1.id)

        tt = gn.distinct_sporozoites_from([(g1,g2)], n_products=[1])
        self.assertEqual(len(tt), 1)
        self.assertEqual(tt[0].parentGenomeIds, (g1.id, g2.id))

        tt = gn.distinct_sporozoites_from([(g1,g2)], n_products=[4])
        self.assertEqual(len(tt), 4)
        self.assertListEqual(np.sum([tx.genome.barcode() for tx in tt], axis=0).tolist(), [2]*gn.num_SNPs())

    def test_distinct(self):
        g1 = gn.Genome.from_reference()
        g2 = gn.Genome.from_barcode([1]*gn.num_SNPs())
        g3 = gn.Genome.from_barcode([1]*6 + [0]*(gn.num_SNPs()-6))
        self.assertListEqual(gn.distinct([g1,g1,g1]), [g1])
        self.assertListEqual(gn.distinct([g1,g2,g3]), [g1,g2,g3])
        self.assertListEqual(gn.distinct([g1,g1,g2]), [g1,g2])

    def test_display(self):
        g3 = gn.Genome.from_barcode([1]*6 + [0]*(gn.num_SNPs()-6))
        self.assertEqual(g3.display_barcode(), '*'*6 + '-'*(gn.num_SNPs()-6))

class TestCrossover(unittest.TestCase):

    nSNP = 24

    def setUp(self):
        self.c1 = np.ones(self.nSNP, dtype=np.uint8)
        self.c2 = np.zeros(self.nSNP, dtype=np.uint8)

    def test_cross(self):
        xpoints = [18]
        c3, c4 = gn.crossover(self.c1, self.c2, xpoints = xpoints)
        print(''.join(map(str,c3)) + '\n' + ''.join(map(str,c4)))
        self.assertListEqual(c3.tolist(), [1]*xpoints[0] + [0]*(self.nSNP-xpoints[0]))
        self.assertListEqual(c4.tolist(), [0]*xpoints[0] + [1]*(self.nSNP-xpoints[0]))

    def test_double_cross(self):
        xpoints = [6, 18]
        c3, c4 = gn.crossover(self.c1, self.c2, xpoints = xpoints)
        print(''.join(map(str,c3)) + '\n' + ''.join(map(str,c4)))
        self.assertListEqual(c3.tolist(), [1]*xpoints[0] + [0]*(xpoints[1]-xpoints[0]) + [1]*(self.nSNP-xpoints[1]))
        self.assertListEqual(c4.tolist(), [0]*xpoints[0] + [1]*(xpoints[1]-xpoints[0]) + [0]*(self.nSNP-xpoints[1]))

if __name__ == '__main__':
    unittest.main()
