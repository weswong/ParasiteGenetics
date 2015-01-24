from Parasite import Parasite
from ParasitePopulation import ParasitePopulation, indices_from_rates

population = ParasitePopulation(n_humans=20, n_parasites=10, n_unique=2)
population.parasites.append(Parasite.from_random())
p=Parasite(123, mixed_mask=255<<16)
population.parasites.append(p)
population.parasites.extend([p.clone() for k in range(5)])
q=population.parasites[-1]
population.parasites.extend([q.outcross(q) for k in range(5)])

print(population)

print('\nBitwise distance, pair counts:\n' + str(population.pairwise_bit_distances()))

print('\nSampling parasite barcodes:\n%s' % population.sample_barcodes(4))

rates = {'clone':0.1, 'outcross':0.05, 'expire':0.2, 'import': 0.5}
import_rate = rates.pop('import',0)
print('\nImport rate=%f' % import_rate)
print('Transitions from rates:\nrates=%s\nindices=%s\n' % (rates, indices_from_rates(rates=rates, n_parasites=len(population.parasites))))

rates = {'clone':0.1, 'outcross':0.05, 'expire':0.2, 'import': 0.5}
print(population.barcode_counts())
population.update(rates)
print(population.barcode_counts())