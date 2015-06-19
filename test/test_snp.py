import math
import unittest

import numpy as np

import genepi.genome as gn
from genepi.snp.snp import SNP
from genepi.event.drug import resistant_sites

class TestBarcode(unittest.TestCase):

    nSNP = 24 # number of barcode SNPs
    nRandom = 100 # number of random allele-frequency tests
    nSigma = 4 # probabilistic cut on binomial allele frequency test

    def setUp(self):
        gn.initialize_from('barcode')

    def test_init(self):
        self.assertEqual(len(gn.Genome.SNP_bins), self.nSNP)
        self.assertEqual(len(gn.Genome.SNP_names), self.nSNP)
        self.assertEqual(len(gn.Genome.SNP_freqs), self.nSNP)
        self.assertEqual(len(gn.Genome.bin_fitness), gn.genome_length())

    def test_reference(self):
        reference = gn.Genome.from_reference()
        self.assertEqual(len(reference.genome), gn.genome_length())
        self.assertListEqual(list(reference.barcode()), [0]*self.nSNP)

    def test_random(self):
        barcode = np.zeros(self.nSNP)
        genomes = [gn.Genome.from_allele_freq() for _ in range(self.nRandom)]
        for i,genome in enumerate(genomes):
            self.assertEqual(len(genome.genome), gn.genome_length())
            barcode += genome.barcode()
        barcode /= self.nRandom
        for f1,f2 in zip(list(barcode),gn.Genome.SNP_freqs):
            self.assertAlmostEqual(f1,f2,delta=self.nSigma*math.sqrt(f2*(1-f2)/self.nRandom))

class TestFullSequence(TestBarcode):

    nSNP = 9311 # number of full-sequence SNPs filtered at 3% and discretized to 1000bp
    nSigma = 5 # increase to 5-sigma for many more allele draws

    def setUp(self):
        gn.initialize_from('sequence', 1000, 0.03)

    def test_add_locus(self):
        chrom, pos = resistant_sites['AM']
        old_snp_bins = gn.Genome.SNP_bins[:]
        add_bin = gn.snp_bin_from_chrom_pos(chrom, pos)
        # SNP already exists at AM-resistance locus after full-genome initialization
        gn.add_locus(chrom, pos, 'Pf.%d.%d' % (chrom, pos), fitness = 0.95)
        self.assertListEqual(gn.Genome.SNP_bins, old_snp_bins)
        self.assertEqual(gn.Genome.bin_fitness[add_bin], 0.95)
        # Test addition of a new SNP to follow a few bins away from    
        offset = 2000
        another_bin = gn.snp_bin_from_chrom_pos(chrom, pos-offset)
        gn.add_locus(chrom, pos-offset, 'Pf.%d.%d' % (chrom, pos-offset), fitness = 0.5)
        self.assertEqual(gn.Genome.SNP_bins[-1], another_bin)
        self.assertEqual(gn.Genome.bin_fitness[another_bin], 0.5)

class TestChromosomes(TestBarcode):

    nSNP = 12 # one per chromosome (none on chromosome #3 or #12)
    bin_size = 1e7 # bigger than longest chromosome

    def setUp(self):
        gn.initialize_from('barcode', self.bin_size)

    def test_chrom_breaks(self):
        self.assertEqual(gn.n_bins_from(5e5), 1)
        self.assertEqual(gn.Genome.chrom_breaks, range(len(gn.chrom_names) + 1))

    def test_bins(self):
        self.assertListEqual(gn.Genome.SNP_bins, [0,1,3,4,5,6,7,8,9,10,12,13])
        self.assertEqual(gn.snp_bin_from_chrom_pos(1,100), 0)
        self.assertEqual(gn.snp_bin_from_chrom_pos(3,100), 2)

def test_closest_snp():
    snps = [SNP(11,5), SNP(11,7), SNP(12,1)]
    assert(gn.find_closest_SNPs(snps) == 2)

def test_rounded_binning():
    assert(gn.get_rounded_binning(12345, ndigits=1) == 10000)
    assert(gn.get_rounded_binning(12345, ndigits=2) == 12000)
    assert(gn.get_rounded_binning(12345, ndigits=3) == 12300)
    assert(gn.get_rounded_binning(12345, ndigits=4) == 12340)

if __name__ == '__main__':
    unittest.main()
