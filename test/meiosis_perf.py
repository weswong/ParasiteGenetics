import itertools
import random
from collections import defaultdict

import numpy as np

from genepi import utils
import genepi.genome as gn

def crossover_list(c1,c2,xpoints):
    c3,c4=c1[:],c2[:]
    if not xpoints:
        return c3,c4
    for l1,l2 in utils.by_pair(xpoints):
        c3[l1:l2], c4[l1:l2] = c4[l1:l2], c3[l1:l2]
    return c3,c4

def crossover_array(c1,c2,xpoints):
    c3=np.copy(c1)
    c4=np.copy(c2)
    if not xpoints:
        return c3,c4
    for l1,l2 in utils.by_pair(xpoints):
        t = np.copy(c3[l1:l2])
        c3[l1:l2] = c4[l1:l2]
        c4[l1:l2] = t
    return c3,c4

def random_list(length):
    return [1 if random.random()<0.5 else 0 for _ in range(length)]

def zeros_list(length):
    return [0 for _ in range(length)]

def ones_list(length):
    return [1 for _ in range(length)]

def zeros_array(length):
    return np.zeros(length,dtype=np.uint8)

def ones_array(length):
    return np.ones(length,dtype=np.uint8)

def random_array(length):
    return np.array(np.random.randint(2, size=length), dtype=np.uint8)

def random_crosspoints(npoints,length):
    return [random.randrange(length) for _ in range(npoints)]

def crossover_test(init1_fn, init2_fn, cross_fn, ntests=1):
    length=100
    c1=init1_fn(length)
    c2=init2_fn(length)
    for _ in range(ntests):
        xpoints=random_crosspoints(5,length)
        c3,c4=cross_fn(c1,c2,xpoints)
    #print(c1,c2,c3,c4)

def all_tests():
    #crossover_test(random_list,random_list,crossover_list)
    #crossover_test(random_array,random_array,crossover_array)
    #crossover_test(zeros_list,ones_list,crossover_list)
    crossover_test(zeros_array,ones_array,crossover_array)

if __name__ == '__main__':
    all_tests()
