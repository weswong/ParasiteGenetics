import warnings
import random
import bisect
import numpy as np # for cumsum, random.poisson
import itertools
from collections import Counter
from Parasite import Parasite

def indices_from_rates(rates, n_parasites):
    indices = { k:[] for k,v in rates.items() }
    rate_names,rate_values = zip(*rates.items())
    n_samples = np.random.poisson(n_parasites * sum(rate_values))
    if n_samples > n_parasites:
        warnings.warn('Requesting more samples (%d) than population (%d).  Defaulting to all parasites.' % (n_samples, n_parasites), RuntimeWarning)
        n_samples=n_parasites
    sampled_idxs = random.sample(range(n_parasites), n_samples)
    rate_cdf = np.cumsum(rate_values)
    rate_cdf /= rate_cdf[-1]
    for i in sampled_idxs:
        k=rate_names[bisect.bisect(rate_cdf,random.random())]
        indices[k].append(i)
    return indices

class ParasitePopulation():
    ''' A collection of clonally unique Parasite objects'''

    def __init__(self, n_humans, n_parasites, propagation_func=lambda:1, n_unique=0):

        if not 0 < n_unique <= n_parasites:
            warnings.warn('Unique barcodes out of range.  Defaulting to every barcode being unique.', RuntimeWarning)
            n_unique = n_parasites

        init_unique_barcodes = [Parasite.random_barcode() for k in range(n_unique)]
        self.parasites = [Parasite(random.choice(init_unique_barcodes)) for k in range(n_parasites)]
        self.n_humans = n_humans

        self.propagation_function = propagation_func

    def __str__(self):
        return 'ParasitePopulation: [\n' + '\n'.join([str(p) for p in self.parasites]) + '\n]'

    def update(self, rates):
        import_rate = rates.pop('import', 0)
        indices = indices_from_rates(rates, len(self.parasites))

        # cache original length for outcrossing selections
        n_parasites = len(self.parasites)

        # append clonal propagation first (to preserve indices for killing)
        for clone_idx in indices.get('clone',[]):
            # draw from overdispersed successful propagation number distribution
            n_repeats = self.propagation_function()
            self.parasites.extend([self.parasites[clone_idx].clone() for k in range(n_repeats)])

        # append outcrossing parasite next
        for outcross_idx in indices.get('outcross',[]):
            # draw from overdispersed successful propagation number distribution
            n_repeats = self.propagation_function()
            self.parasites.extend( [ self.parasites[outcross_idx].outcross( self.parasites[ random.randrange(n_parasites) ] ) for k in range(n_repeats) ] )

        # append import cases
        if import_rate:
            n_imports = np.random.poisson(import_rate)
            for import_idx in range(n_imports):
                self.parasites.append(Parasite.from_random())

        # expire parasites last, removing from back to front to preserve indices
        for expire_idx in sorted(indices.get('expire',[]), reverse=True):
            # TODO: performance?
            del self.parasites[expire_idx]

    def barcodes(self):
        return [p.barcode for p in self.parasites]

    def infections(self):
        return [(p.barcode,p.mixed_mask) for p in self.parasites]

    def mean_COI(self):
        return len(self.parasites) / float(self.n_humans)

    def sample_barcodes(self, N):
        n_parasites=len(self.parasites)
        if N > n_parasites:
            warnings.warn('Requesting sample (%d) larger than population (%d).  Defaulting to all parasites.' % (N, n_parasites), RuntimeWarning)
            N = n_parasites
        return [p.barcode for p in random.sample(self.parasites, N)]

    def barcode_counts(self):
        counts = Counter()
        for p in self.parasites:
            counts[p.barcode] += 1
        return counts

    def pairwise_bit_distances(self):
        bit_distances = Counter()
        for p1,p2 in itertools.combinations(self.parasites,2):
            bit_distances[Parasite.bitwise_distance(p1.barcode,p2.barcode)] += 1
        return bit_distances