import itertools
import random
from collections import defaultdict
import copy
import numpy as np
import hashlib


def zeros_array(length):
    return np.zeros(length,dtype=np.uint8)

def ones_array(length):
    return np.ones(length,dtype=np.uint8)

def random_array(length):
    return np.array(np.random.randint(2, size=length), dtype=np.uint8)

#@profile
def genome_hash(genome):
    V=(tuple(v) for v in genome.values())
    T=tuple(V)
    return hash(T)

@profile
def genome_hash2(genome):
    T=tuple(genome)
    return hash(T)

@profile
def genome_hash3(genome):
    v = genome.view(np.uint8)
    return hashlib.sha1(v)#.hexdigest()

def distinct(genomes,hash_fn):
    distinct=[]
    seen = set()
    for g in genomes:
        h=hash_fn(g)
        if h not in seen:
            seen.add(h)
            distinct.append(g)
    return distinct

def distinct_test(hash_fn,ntests=1):
    length=100
    genomes=[]
    genomes.append({1:random_array(length),2:random_array(length)})
    genomes.append(copy.deepcopy(genomes[0]))
    genomes.append({1:random_array(length),2:random_array(length)})
    for _ in range(ntests):
        D=distinct(genomes,hash_fn)
        #print('%d distinct genomes out of %d'%(len(D),len(genomes)))

@profile
def distinct_test2(hash_fn,ntests=1):
    length=100
    genomes=[]
    genomes.extend([random_array(length) for _ in range(3)])
    genomes.append(copy.deepcopy(genomes[0]))
    genomes.append(copy.deepcopy(genomes[2]))
    for _ in range(ntests):
        D=distinct(genomes,hash_fn)
        #print('%d distinct genomes out of %d'%(len(D),len(genomes)))

if __name__ == '__main__':
    #distinct_test(genome_hash)
    #distinct_test2(genome_hash2,100)
    distinct_test2(genome_hash3,100)
