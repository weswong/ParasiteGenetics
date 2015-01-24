import os
import json
import math
import collections
import pprint
from scipy import stats
import warnings
import numpy as np
from calculate_sim_metrics import BarcodeMetricCalculator

class BarcodeLikelihoodComparator():

    data_unique_counts = [ 84, 109, 59, 57,  46,  55,  89, 161 ]
    data_mono_samples  = [ 90, 138, 96, 81, 100, 112, 158, 236 ]
    data_poly_samples  = [ 41,  68, 13, 15,  26,  26, 110,  78 ]

    data_repeats_equal_two   = [ 3, 9, 6, 6, 9, 6,  7,  9 ]
    data_repeats_greater_two = [ 0, 3, 4, 3, 6, 6, 10, 12 ]

    persist_years = [ 21, 5, 3, 1 ]
    created_persists = 28
    destroyed_persists = 18

    data = { 'mixed_infection_fraction' : [ p/float(p+m) for m,p in zip(data_mono_samples,data_poly_samples) ],
             'unique_fraction'          : [ u/float(m) for u,m in zip(data_unique_counts,data_mono_samples) ],
             'persist_years'            : persist_years,
             'created_persists'         : created_persists,
             'destroyed_persists'       : destroyed_persists,
             'repeats_equal_two'        : data_repeats_equal_two,
             'greater_two_per_repeat'   : [ g/float(e+g) for e,g in zip(data_repeats_equal_two,data_repeats_greater_two) ] }

    data_err = { 'mixed_infection_fraction' : [ math.sqrt(f*(1-f)/(m+p)) for f,m,p in zip(data['mixed_infection_fraction'],data_mono_samples,data_poly_samples) ],
                 'unique_fraction'          : [ math.sqrt(f*(1-f)/m) for f,m in zip(data['unique_fraction'],data_mono_samples) ],
                 'persist_years'            : [ math.sqrt(p) for p in data['persist_years'] ],
                 'created_persists'         : math.sqrt(data['created_persists']),
                 'destroyed_persists'       : math.sqrt(data['destroyed_persists']),
                 'repeats_equal_two'        : [ math.sqrt(r) for r in data['repeats_equal_two'] ],
                 'greater_two_per_repeat'   : [ math.sqrt(f*(1-f)/(e+g)) for f,e,g in zip(data['greater_two_per_repeat'],data_repeats_equal_two,data_repeats_greater_two) ] }

    # correct p(1-p) error approximation for one bin where p=0
    data_err['greater_two_per_repeat'][0] = 0.24788 # 68% of beta(a=1,b=4) below this value

    ###############################################
    # SYSTEMATIC EVALUATION OF DATA-ANALYSIS CUTS #
    ###############################################

    # Changing from N>1 to N>0 as polygenomic definition
    strict_polygenomic = {
    'unique_counts' : [ 81,  98, 57, 56, 43,  48,  66, 129 ],
    'mono_samples'  : [ 87, 130, 94, 80, 95, 108, 135, 198 ],
    'poly_samples'  : [ 44,  76, 15, 16, 31,  30, 133, 116 ]
    }

    # Loosening haplotyping from counting N and X as mismatch
    # to allowing at least one of each while still calling the same
    loose_haplotyping = {
    'unique_counts' : [ 81,  98, 59, 55,  41,  46,  78, 126 ],
    'mono_samples'  : [ 90, 138, 96, 81, 100, 112, 158, 236 ],
    'repeats_equal_two'   : [ 3, 11, 6, 7, 7, 5,  7, 12 ],
    'repeats_greater_two' : [ 1,  4, 3, 3, 7, 8, 12, 16 ],
    'persist_years' : [ 27, 6, 7, 1 ],
    'created_persists' : 38,
    'destroyed_persists' : 26
    }

    mixed_infection_fraction_strict = [ p/float(p+m) for m,p in zip(strict_polygenomic['mono_samples'],strict_polygenomic['poly_samples']) ]
    unique_fraction_strict          = [ u/float(m) for u,m in zip(strict_polygenomic['unique_counts'],strict_polygenomic['mono_samples']) ]
    unique_fraction_loose           = [ u/float(m) for u,m in zip(loose_haplotyping['unique_counts'],loose_haplotyping['mono_samples']) ]
    greater_two_loose               = [ g/float(e+g) for e,g in zip(loose_haplotyping['repeats_equal_two'],loose_haplotyping['repeats_greater_two']) ]

    def combined_means_and_errors(d,de,metric, alt_metric_list):
        #print(metric,'Before',d[metric],de[metric])
        # Take mean of different cuts for value
        # Add in quadrature half the min-max spread to the statistical uncertainty
        sys_err = 0.5 * ( np.nanmax( [d[metric]]+alt_metric_list, axis=0 )-np.nanmin( [d[metric]]+alt_metric_list, axis=0 ) )
        #print(metric,'sys_err',sys_err)
        if isinstance(d[metric], collections.Iterable):
            d[metric]  = [ sum(D)/(1.0+len(alt_metric_list)) for D in zip(d[metric],*alt_metric_list) ]
            de[metric] = [ math.sqrt(a**2+b**2) for (a,b) in zip( de[metric], sys_err) ]
        else:
            d[metric]  = sum([d[metric]]+alt_metric_list)/(1.0+len(alt_metric_list))
            de[metric] = math.sqrt(de[metric]**2+sys_err**2)
        #print(metric,'After',d[metric],de[metric])

    combined_means_and_errors(data, data_err, 'mixed_infection_fraction', [mixed_infection_fraction_strict])
    combined_means_and_errors(data, data_err, 'unique_fraction', [unique_fraction_strict,unique_fraction_loose])
    combined_means_and_errors(data, data_err, 'persist_years', [loose_haplotyping['persist_years']])
    combined_means_and_errors(data, data_err, 'created_persists', [loose_haplotyping['created_persists']])
    combined_means_and_errors(data, data_err, 'destroyed_persists', [loose_haplotyping['destroyed_persists']])
    combined_means_and_errors(data, data_err, 'repeats_equal_two', [loose_haplotyping['repeats_equal_two']])
    combined_means_and_errors(data, data_err, 'greater_two_per_repeat', [greater_two_loose])

    def __init__(self, work_dir='simulations', idxs=[0]):
        self.calculators=[BarcodeMetricCalculator(work_dir,i) for i in idxs]
        self.compare(verbose=False)#True)

    def __str__(self):
        pp = pprint.PrettyPrinter(indent=4)
        return pp.pformat(self.chi_squares)

    def compare(self, verbose=False):

        self.chi_squares={}

        def chi2(x):
            k,d,s,de,se = x
            xs=(d-s)**2 / (de**2+se**2)
            if verbose: print(d,s,de,se,xs)
            self.chi_squares[k].append(xs)

        warnings.simplefilter('once', UserWarning)
        for k,v in self.data.items():
            if verbose: print(k)
            self.chi_squares[k]=[]
            sim_errs = [calc.sim_err.get(k,[]) for calc in self.calculators]
            sim_data = [calc.sim.get(k,[]) for calc in self.calculators]
            sim_mean = np.mean(sim_data,0)

            if len(sim_data)<=1 and not any(sim_errs):
                warnings.warn('Using data uncertainty for sim variation is likely an underestimate.')
                use_sim_err=False
                #print(sim_errs,sim_data,sim_mean)

            else:
                use_sim_err=True
                sim_std  = np.std(sim_data,0) if len(sim_data) > 1 else 0
                sim_errs_mean = np.mean(sim_errs,0) if sim_errs[0] else 0
                sim_err = np.sqrt(sim_errs_mean**2+sim_std**2)
                #print(sim_errs,sim_data,sim_mean,sim_std,sim_errs_mean,sim_err)

            if isinstance(v, collections.Iterable):
                for i,vi in enumerate(v):
                    chi2((k,vi,sim_mean[i],self.data_err[k][i],sim_err[i] if use_sim_err else self.data_err[k][i]))
            else:
                chi2((k,v,sim_mean,self.data_err[k],sim_err if use_sim_err else self.data_err[k]))

    def get_likelihood(self):
        chi2,ndof=0,0
        for v in self.chi_squares.values():
            #ndof += len(v)
            #chi2 += sum(v)
            ndof += 1
            chi2 += sum(v)/len(v)
        prob= 1.0-stats.chi2.cdf(chi2, ndof)
        return {'chi2':chi2,'ndof':ndof,'prob':prob,'LL':-0.5*chi2}

if __name__ == '__main__':
    work_dir_full=os.path.join('simulations/test_likelihoods', 'sim_%d'%8)
    nseeds=10
    
    #for i in range(nseeds):
    #    cmp=BarcodeLikelihoodComparator(work_dir=work_dir_full,idxs=range(i,i+1))
    #    print(cmp)
    #    print(cmp.get_likelihood())

    cmp=BarcodeLikelihoodComparator(work_dir=work_dir_full,idxs=range(nseeds))
    print(cmp)
    print(cmp.get_likelihood())
    with open(os.path.join(work_dir_full,'Likelihoods.json'),'w') as likelihood_file:
        json.dump({'summary':cmp.get_likelihood(),'chi_squares':cmp.chi_squares},likelihood_file)


