import os
import json
import sys
from ParasiteSimulation import ParasiteSimulation
from functools import reduce
from compare_likelihood import BarcodeLikelihoodComparator

def run_simulation(work_dir, nseeds=1, param_override_fn=lambda x:None):
    dirname=os.path.dirname(os.path.abspath(__file__))
    cfg_path = os.path.join(dirname,'config.json')
    cp_path = os.path.join(work_dir, 'config.json') # TODO: generalize cfg name

    with open(cfg_path,'r') as cfg_f:
        cfg=json.loads(cfg_f.read())

    param_override_fn(cfg)

    if not os.path.exists(os.path.dirname(cp_path)):
        os.makedirs(os.path.dirname(cp_path))

    with open(cp_path,'w') as cp_f:
        json.dump(cfg, cp_f)

    for seed in range(nseeds):
        sim = ParasiteSimulation.from_config_file(cp_path, seed=seed)
        sim.working_directory = work_dir

        print(sim.params)
        sim.run()

def getFromDict(dataDict, mapList):
    return reduce(dict.__getitem__, mapList, dataDict)

def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value

def override(param_overrides):
    def f(cfg):
        for p in param_overrides:
            name,value = p.split(':')
            #cfg[name] = float(value)
            setInDict(cfg, name.split('.'), float(value)) # TODO: generalize to non-float values and integer list indices
    return f

if __name__ == "__main__":

    nargs=len(sys.argv)

    if nargs == 1:
        run_simulation(work_dir='simulations')
    else:
        work_dir = sys.argv[1]
        idx = sys.argv[2]
        print('Index: %s' % idx)
        nseeds = int(sys.argv[3])
        print('# of random seeds: %d' % nseeds)
        params = sys.argv[4:]
        print('Command line params: %s' % params)
        #workdir_description = '_'.join([p.replace(':','_') for p in params])
        workdir_description = 'sim_' + idx
        
        work_dir_full = os.path.join(work_dir, workdir_description)

        # TODO: remove the block below, which is now redundant with the sim_<idx> naming above?
        # Write index.txt after first making directory
        if not os.path.exists( work_dir_full ):
            os.makedirs( work_dir_full )
        
        index_fn = os.path.join(work_dir_full, 'CalibPointsIndex.txt')
        print(index_fn)
        with open( index_fn, 'w' ) as index_file:
            index_file.write(idx)
        
        run_simulation( work_dir = work_dir_full,
                        nseeds = nseeds,
                        param_override_fn=override(params) )

        cmp=BarcodeLikelihoodComparator(work_dir=work_dir_full,idxs=range(nseeds))
        print(cmp)
        print(cmp.get_likelihood())
        with open(os.path.join(work_dir_full,'Likelihoods.json'),'w') as likelihood_file:
            json.dump({'summary':cmp.get_likelihood(),'chi_squares':cmp.chi_squares},likelihood_file)
