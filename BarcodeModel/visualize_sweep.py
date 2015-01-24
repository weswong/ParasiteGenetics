import os
import json
import glob

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

from run_simulation import getFromDict

def visualize(work_dir,param_names):
    point_chi2s=[]
    for dir in glob.glob('%s/sim_*' % work_dir):
        with open(os.path.join(dir,'config.json')) as cf:
            params=json.loads(cf.read())
        param_values=[getFromDict(params,n.split('.')) for n in param_names]
        with open(os.path.join(dir,'likelihoods.json')) as lf:
            chi2=json.loads(lf.read())['summary']['chi2']
        #print(param_names,param_values,chi2)
        point_chi2s.append(param_values+[chi2])

    points_and_chi2s=zip(*point_chi2s)
    nparams=len(param_names)

    if nparams==1:
        plt.scatter(*points_and_chi2s)

    elif nparams==2:
        x = np.array(points_and_chi2s[0])
        y = np.array(points_and_chi2s[1])
        z = np.array(points_and_chi2s[2])
        xi,yi = np.mgrid[x.min():x.max():500j,y.min():y.max():500j]

        def normal_interp(x, y, a, xi, yi):
            rbf = interpolate.Rbf(x, y, a, 
                                  function='cubic')
            ai = rbf(xi, yi)
            return ai

        def rescaled_interp(x, y, a, xi, yi):
            a_rescaled = (a - a.min()) / a.ptp()
            ai = normal_interp(x, y, a_rescaled, xi, yi)
            ai = a.ptp() * ai + a.min()
            return ai

        def plot(x, y, a, ai, title):
            fig, ax = plt.subplots(num='chi2_interpolation')
            cmap='jet_r'
            vmin,vmax=5,20
            im = ax.imshow(ai.T, origin='lower',
                           extent=[x.min(), x.max(), y.min(), y.max()],
                           aspect=(x.max()-x.min())/(y.max()-y.min()),
                           vmin=vmin, vmax=vmax, cmap=cmap)
            ax.scatter(x, y, c=a, 
                       vmin=vmin, vmax=vmax, cmap=cmap)

            ax.set(xlabel=param_names[0], 
                   ylabel=param_names[1], 
                   title=title)
            fig.colorbar(im)

        z_rescale = rescaled_interp(x, y, z, xi, yi)
        plot(x, y, z, z_rescale,r'$\chi^2$ interpolation')

    else:
        raise Exception("Haven't implemented 3-d plotting yet")

    plt.show()

if __name__ == '__main__':
    visualize('simulations/R0_rebound_year_R0_rebound_sweep',
              ['R0_rebound_year','R0_rebound'])