import os
import time
import threading
import itertools
import json
import numpy as np
import subprocess
import sys

class ThreadManager(threading.Thread):
    lock = threading.Lock()
    num_threads = 0

    def __init__(self, args, work_dir, idx, start_time, nseeds=1):
        threading.Thread.__init__(self)
        self.Args = args
        self.Idx = idx
        self.nSeeds = nseeds
        self.Start_Time = start_time
        self.WorkDir = work_dir

    def run(self):

        ThreadManager.lock.acquire()      # Begin critical section
        ThreadManager.num_threads += 1
        ThreadManager.lock.release()      # End critical section

        dirname=os.path.dirname(os.path.abspath(__file__))
        subprocess.call(['python',os.path.join(dirname,'run_simulation.py')] + [self.WorkDir] + [str(self.Idx)] + [str(self.nSeeds)] + self.Args)

        ThreadManager.lock.acquire()      # Begin critical section
        ThreadManager.num_threads -= 1
        ThreadManager.lock.release()      # End critical section

        print("Ending " + str(self.Idx) + " for " + str(self.Args) + ", Elapsed time: " + str( time.time() - self.Start_Time ))

def parameter_sweep(param_names, param_combos, work_dir='.', nseeds=3, num_cpus=4):

    threads = []
    start_time = time.time()

    for (idx, params) in enumerate(param_combos):

        while ThreadManager.num_threads >= num_cpus:
            time.sleep(1)

        args = [':'.join([n,str(v)]) for (n,v) in zip(param_names,params)] # TODO: enforce floating-point string format?

        print("Starting " + str(idx) + " for params: " + str(args))

        thread = ThreadManager(args, work_dir, idx, time.time(), nseeds)
        threads += [thread]
        thread.start()
        time.sleep(1)

    for t in threads:
        t.join()

    print("Elapsed time: " + str( time.time() - start_time ))

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print('No path to parameter_values specified.  Defaulting to sweep.')

        #param_ranges = {
        #    'R0_initial' : np.arange(1.0, 5.0, 0.2),
        #    'R0_reduction_scale' : np.arange(0.2, 1.4, 0.1) }

        param_ranges = {
            'R0_initial' : [1.9, 2.1],
            'R0_final' : [1.4, 1.5],
            'rates.expire' : [0.02, 0.03] }

        param_names = param_ranges.keys()
        param_combos = itertools.product(*param_ranges.values())
        work_dir = 'simulations'

    else:
        if not os.path.exists(sys.argv[1]):
            sys.exit('ERROR: File %s was not found!' % sys.argv[1])

        # load parameter combos from file
        tmpjson = sys.argv[1]
        with open(tmpjson,'r') as f:
            j=json.loads(f.read())
        param_names = j['param_names']
        param_combos = j['param_combos']
        work_dir = j['work_dir']

    parameter_sweep(param_names, param_combos, work_dir=work_dir)
