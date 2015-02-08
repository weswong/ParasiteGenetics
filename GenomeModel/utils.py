import logging
logging.basicConfig(format='%(message)s', level=logging.DEBUG)
log = logging.getLogger(__name__)

def choose_without_replacement(from_list,N):
    '''
    O(M) in choose M from N scenario,
    which is much faster for typical use case
    than random.choice, which is O(N)
    '''
    chosen_idxs=set()
    chosen=[]
    size=len(from_list)
    for j in range(size-N,size):
        t = random.randint(j) # 0 to j-1
        idx = t if t not in chosen_idxs else j
        chosen_idxs.add(idx)
        chosen.append(from_list[idx])
    return chosen
