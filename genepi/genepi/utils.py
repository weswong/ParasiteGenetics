import math
import random
import itertools

def pairwise(l):
    "l -> (l0,l1), (l2,l3), (l4, l5), ..."
    if len(l) % 2:
        l.append(None)
    a = iter(l)
    return itertools.izip(a, a)

def accumulate_cdf(iterable):
    cdf,subtotal=[],0
    norm=float(sum(iterable))
    if not norm:
        return []
    for it in iterable:
        subtotal += it/norm
        cdf.append(subtotal)
    cdf[-1]=1.0
    return cdf

def nextTime(rateParameter):
    if rateParameter<=0:
        raise Exception('Cannot calculate next time from zero rate.')
    return -math.log(1.0 - random.random()) / rateParameter

def poissonRandom(rateParameter):
    if rateParameter<=0:
        return 0
    sumTime=0
    N=0
    while True:
        sumTime+=nextTime(rateParameter)
        if sumTime>1:
            break
        N+=1
    return N

def weighted_choice(cumwts):
    R = random.random()
    idx=sum(itertools.takewhile(bool, (cw < R for cw in cumwts)))
    return idx

def choose_with_replacement(M,N):
    choose = random.choice
    indices = range(N)
    return [choose(indices) for _ in range(M)]

def choose_without_replacement(M,N):
    '''
    Floyd algorithm: O(M) in choose M from N scenario,
    which can be much faster for some typical use cases
    than random.sample, which seems to be O(N)
    '''
    if M>N:
        raise Exception('Cannot sample %d from %d without replacement',(M,N))
    if M==N:
        return range(M)
    chosen_idxs=set()
    for j in range(N-M,N):
        t = random.randint(0,j)
        idx = t if t not in chosen_idxs else j
        chosen_idxs.add(idx)
    return list(chosen_idxs)
