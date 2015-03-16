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

#@profile
def genome_hash2(genome):
    T=tuple(genome)
    return hash(T)

#@profile
def genome_hash3(genome):
    v = genome.view(np.uint8)
    return hashlib.sha1(v)#.hexdigest()

@profile
def distinct(genomes,hash_fn):
    distinct=[]
    seen = set()
    for g in genomes:
        h=hash_fn(g)
        if h not in seen:
            seen.add(h)
            distinct.append(g)
    return distinct

@profile
def distinct2(genomes):
    b = np.ascontiguousarray(genomes).view(np.dtype((np.void, genomes.dtype.itemsize * genomes.shape[1])))
    #_, idx = np.unique(b, return_index=True)
    #unique_genomes = genomes[idx]
    unique_genomes = np.unique(b).view(genomes.dtype).reshape(-1, genomes.shape[1])
    return unique_genomes
    #return np.unique(np.ascontiguousarray(genomes).view(np.dtype((np.void, genomes.dtype.itemsize * genomes.shape[1])))).view(genomes.dtype).reshape(-1, genomes.shape[1])

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
    length=100000
    genomes=[]
    genomes.extend([random_array(length) for _ in range(3)])
    genomes.append(copy.deepcopy(genomes[0]))
    genomes.append(copy.deepcopy(genomes[2]))
    for _ in range(ntests):
        D=distinct(genomes,hash_fn)
        #print('%d distinct genomes out of %d'%(len(D),len(genomes)))

@profile
def distinct_test3(ntests=1):
    length=100000
    genomes=[]
    genomes.extend([random_array(length) for _ in range(3)])
    genomes.append(copy.deepcopy(genomes[0]))
    genomes.append(copy.deepcopy(genomes[2]))
    genomes=np.asarray(genomes)
    for _ in range(ntests):
        #genomes=np.asarray(genomes)
        #print(genomes)
        D=distinct2(genomes)
        #print(D)

if __name__ == '__main__':
    #distinct_test(genome_hash)
    #distinct_test2(genome_hash2,100)
    distinct_test2(genome_hash3,100)
    distinct_test3(100)
