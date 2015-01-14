import json
import random

class Parasite():
    ''' A clonally unique infection holding a parasite genetic barcode'''
    numSites = 24

    def __init__(self, barcode, mixed_mask=0):
        self.barcode=barcode
        self.mixed_mask=mixed_mask

    def __str__(self):
        s = 'Parasite: ' + ('{0:0>%db}' % self.numSites).format(self.barcode)
        if self.mixed_mask:
            s += ', %d mixed = ' % self.n_mixed_positions() + ('{0:0>%db}' % self.numSites).format(self.mixed_mask)
        return s

    @classmethod
    def random_barcode(cls):
        return random.getrandbits(cls.numSites)

    @classmethod
    def from_random(cls):
        barcode = cls.random_barcode()
        return cls(barcode)

    @classmethod
    def from_outcrossing(cls, barcode1, barcode2):
        mask = cls.random_barcode()
        barcode = (barcode1 & mask) | (barcode2 & ~mask)
        return cls(barcode)

    @staticmethod
    def bitwise_distance(barcode1, barcode2):
        return bin(barcode1 ^ barcode2).count('1')

    def n_mixed_positions(self):
        return bin(self.mixed_mask).count('1')

    def clone(self):

        if not self.mixed_mask:
            return Parasite(self.barcode)

        # Zero-inflated binomial process of having some chance
        # to lose all variability in one generation
        if random.random() < 0.2:
            # Determine which mixed sites will change in new barcode
            change_mask = self.mixed_mask & self.random_barcode()
            return Parasite(self.barcode ^ change_mask)

        # Otherwise, some per-position rate of losing mixed positions
        # TODO: if we want to get more refined than 50% per site
        #       we should change how the following line is done
        mixed_mask = self.mixed_mask & self.random_barcode()
        # Determine which mixed sites will change in new barcode
        # N.B. This will be randomly flipping bits at both positions
        # that lose variability and those that retain it.
        change_mask = self.mixed_mask & self.random_barcode()
        return Parasite(self.barcode ^ change_mask, mixed_mask)

    def outcross(self, p2):
        # Sites where the two barcodes differ or where either one separately has a mixed position
        variable_positions = (self.barcode ^ p2.barcode) | self.mixed_mask | p2.mixed_mask
        return Parasite(self.barcode,variable_positions).clone()