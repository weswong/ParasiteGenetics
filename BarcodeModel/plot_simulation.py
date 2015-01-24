import os
import json
import math
import random
from itertools import chain
from Parasite import Parasite
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from plotting import adjust_spines

def barcode2color(barcode):
    B_MASK = 255
    G_MASK = 255<<8
    R_MASK = 255<<16

    r = (barcode & R_MASK)>>16
    g = (barcode & G_MASK)>>8
    b = barcode & B_MASK

    #return [c/256. for c in [r,g,b]]
    #return [0.5*(0.7+c/256.) for c in [r,g,b]] # mix 50-50 with light neutral color for more pleasing palette?
    return [0.5*(0.9+c/256.) for c in [0.2*r,0.5*g,b]] # mix 50-50 with light neutral color for more pleasing palette?

def plot_simulation(work_dir='simulations', idx=0):

    #work_dir = os.path.join('simulations','sweep_import_0.100000_R0_4.000000')
    print('Plotting from output directory: %s' % work_dir)

    with open(os.path.join(work_dir,'barcode_report_%d.json' % idx),'r') as infile:
        bcj = json.loads(infile.read())

    n_parasites = bcj['n_parasites']
    n_barcodes = bcj['n_barcodes']
    n_mixed = bcj['n_mixed']
    p_outcross = bcj['p_outcross']
    fraction_mixed = [[m/float(p) if p>0 else 0 for (m,p) in zip(my,py)] for (my,py) in zip(n_mixed,n_parasites)]
    census = bcj['barcode_census']
    n_years = len(census)

    fig = plt.figure('Parasites + Barcodes', figsize=(12,6), facecolor='w')
    ax = fig.add_subplot(211)
    for pop_idx in range(len(n_parasites)):
        plt.plot(n_parasites[pop_idx],'midnightblue')
        plt.plot(n_barcodes[pop_idx],'skyblue')
        #plt.plot(n_mixed[pop_idx],'darkseagreen')
        if not pop_idx:
            leg=plt.legend(('# of parasite infections',
                            '# of unique parasite barcodes', 
                            #'# of mixed infections'
                            ), 
                           frameon=False)

    for l in leg.legendHandles:
        l.set_linewidth(10)

    adjust_spines(ax,['left','bottom'])
    ax.set_xticks(range(0, 365*n_years, 2*365))
    ax.set_xticklabels(range(2013-n_years+1, 2013+1, 2))

    ax = fig.add_subplot(212)
    for pop_idx in range(len(n_parasites)):
        plt.plot(p_outcross[pop_idx],'darkolivegreen')
        plt.plot(fraction_mixed[pop_idx],'darkseagreen')
        if not pop_idx:
            leg=plt.legend(('outcrossing probability','mixed-infection fraction'), frameon=False)

    for l in leg.legendHandles:
        l.set_linewidth(10)

    adjust_spines(ax,['left','bottom'])
    ax.set_xticks(range(0, 365*n_years, 2*365))
    ax.set_xticklabels(range(2013-n_years+1, 2013+1, 2))

    plt.tight_layout()
    plt.savefig(os.path.join(work_dir,'parasites_and_barcodes.png'))

    ####################
    annual_samples = [100]*n_years
    annual_samples[-8:] = [90, 138, 96, 81, 100, 112, 158, 236] # using Thies sampled numbers in last 8 years
    ####################

    last_year = 2013

    N = float(Parasite.numSites)

    allele_freqs = []

    repeats = {}
    n_mixed_positions = []
    n_random_pairs_mixed_positions = []

    patch_handles = []
    left = np.zeros(n_years,)
    row_counts = np.zeros(n_years,)
    fig = plt.figure('Repeat barcodes', figsize=(10,8), facecolor='w')
    ax = fig.add_subplot(111)

    def plot_bar(year,value,norm,color,textcolor='w'):
        patch_handles.append(ax.barh(n_years-year, value/norm, align='center', left=left[year], color=color))
        left[year] += value/norm
        row_counts[year] += 1
        # we know there is only one patch but could enumerate if expanded
        patch = patch_handles[-1][0] 
        bl = patch.get_xy()
        x = 0.5*patch.get_width() + bl[0]
        y = 0.5*patch.get_height() + bl[1]
        if value != 2:
            ax.text(x, y, "%d" % value, ha='center',va='center', color=textcolor)

    def allele_frequencies_from_samples(samples):
        freqs=np.zeros(N)
        for s in samples:
            freqs += [int(d) for d in ('{0:0>%db}' % N).format(s[0])]
        freqs /= float(len(samples))
        return freqs

    random.seed(1234)
    for year in range(7,n_years):
        annual_infections = census[str(year)]
        '''
        if len(annual_infections) <= annual_samples[year]:
            sampled_infections = annual_infections
        else:
            sampled_infections = random.sample(annual_infections, annual_samples[year])
        mixed_samples = [s for s in sampled_infections if s[1]]
        unmixed_samples = [s for s in sampled_infections if not s[1]]
        n_mixed_positions.append([bin(s[1]).count('1') for s in mixed_samples])

        allele_freqs.append(allele_frequencies_from_samples(unmixed_samples))
        '''

        # another sampling that chooses N unmixed samples rather than the unmixed set out of N total samples
        all_unmixed = [s for s in annual_infections if not s[1]]
        unmixed_samples = random.sample(all_unmixed, annual_samples[year])

        '''
        rand_diff_bits = []
        for i in range(1000):
            s1=random.choice(unmixed_samples)
            s2=random.choice(unmixed_samples)
            rand_diff_bits.append(bin(s1[0] ^ s2[0]).count('1'))

        n_random_pairs_mixed_positions.append(rand_diff_bits)
        '''

        barcode_counter=Counter([p[0] for p in unmixed_samples])
        norm=float(len(unmixed_samples))
        singles=0
        for k,v in sorted(barcode_counter.items(), key=lambda tup: tup[1]):
            if v==1:
                singles += 1
                continue
            else:
                if singles:
                    plot_bar(year, singles, norm, '#eeeeee', '#888888')
                    singles=0
            #plot_bar(year, v, norm, '#00709a')        # grayish blue for all
            plot_bar(year, v, norm, barcode2color(k))  # color to visualize year-to-year persistence
        if singles:
            plot_bar(year, singles, norm, '#eeeeee', '#888888')

        for k in repeats.keys():
            if k not in barcode_counter:
                repeats[k].append(0)

        for k,v in barcode_counter.items():
            if k not in repeats.keys():
                repeats[k]=[0]*year
            repeats[k].append(v)

    ax.set_yticks(range(1,9))
    ax.set_yticklabels(range(2013,2013-n_years,-1))
    ax.get_xaxis().set_visible(False)
    for t in ax.yaxis.get_ticklines():
        t.set_visible(False)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    plt.xlim([0,1])
    plt.ylim([-0.02,9])
    plt.tight_layout()
    plt.savefig(os.path.join(work_dir,'repeat_barcodes.png'))

    def n_initial_zeros(ll):
        n=0
        for l in ll:
            if l==0:
                n+=1
            else:
                break
        return n

    i,xx,yy,ss,ll,cc=0,[],[],[],[],[]
    fig = plt.figure('Persistent barcodes', figsize=(9,10), facecolor='w')
    ax = fig.add_subplot(111)
    #for k,v in repeats.items():
    for k,v in sorted(repeats.items(), key=lambda t:(n_initial_zeros(t[1]),t[0]), reverse=True):
        if sum([n>0 for n in v]) <= 1: # or sum(v) <= 2:
            continue
        xx.extend(range(n_years))
        yy.extend([i]*n_years)
        ss.extend(v)
        ll.extend([str(k)])
        cc.extend([barcode2color(k)]*n_years)
        nonzero=[t for t, e in enumerate(v) if e != 0]
        ax.plot([nonzero[0],nonzero[-1]],[i]*2,'k-',alpha=0.2)
        i+=1

    #scatter_color = 'gray' # gray for all bubbles
    scatter_color = cc      # color to match within-year repeat bars
    ax.scatter(xx, yy, [40*sz for sz in ss], c=scatter_color, linewidth=0.1)
    plt.ylim([-0.5, len(ll)-0.5])
    plt.xlim([-0.5, n_years-0.5])
    adjust_spines(ax,['left','bottom'])
    ax.set_xticks(range(0, n_years))
    ax.set_xticklabels(range(2013-n_years+1, 2013+1))
    ax.set_yticks(range(len(ll)))
    ax.set_yticklabels(ll)
    ax.set_xlim([n_years-8.1,n_years])
    plt.tight_layout()
    plt.savefig(os.path.join(work_dir,'persistent_barcodes.png'))

    '''
    fig=plt.figure('Mixed-infection similarity',figsize=(14,8),facecolor='w')
    gs1 = gridspec.GridSpec(2,1)
    gs1.update(left=0.06, right=0.48, top=0.96, hspace=0.25)
    ax=plt.subplot(gs1[0,0])
    ax.text(0,0.95,'Mixed infections', transform=ax.transAxes, fontsize=14)
    mixed_similarity = [ 1.-n/N for n in chain(*n_mixed_positions)]
    plt.hist(mixed_similarity, bins=np.arange(-0.5/N,1+0.5/N,1/N), width=0.8/N, color='thistle')
    adjust_spines(ax,['left','bottom'])
    ax.set_xticks(np.arange(0, 1, 0.2))
    plt.xlim([-0.5/N,1+0.5/N])

    ax=plt.subplot(gs1[1,0])
    ax.text(0,0.95,'Random pairs of single infections', transform=ax.transAxes, fontsize=14)
    mixed_similarity = [ 1.-n/N for n in chain(*n_random_pairs_mixed_positions)]
    plt.hist(mixed_similarity, bins=np.arange(-0.5/N,1+0.5/N,1/N), width=0.8/N, color='lightskyblue')
    adjust_spines(ax,['left','bottom'])
    ax.set_xticks(np.arange(0, 1, 0.2))
    plt.xlim([-0.5/N,1+0.5/N])
    plt.xlabel('Fraction of shared bits', fontsize=12)

    gs2 = gridspec.GridSpec(4,2)
    gs2.update(left=0.55, right=0.98, top=0.96, bottom=0.06, wspace=0.2, hspace=0.2)

    for i,year in enumerate(range(-8,0,2)):
        ax=plt.subplot(gs2[i])
        ax.text(0,0.8,str(2013+year+1) + ' - ' + str(2013+year+2), fontsize=11, transform=ax.transAxes)
        mixed_similarity = [ 1.-n/N for n in n_mixed_positions[year]]
        plt.hist(mixed_similarity, bins=np.arange(-0.5/N,1+0.5/N,1/N), width=0.8/N, color='thistle')
        adjust_spines(ax,[])
        plt.xlim([-0.5/N,1+0.5/N])

    for i,year in enumerate(range(-8,0,2)):
        ax=plt.subplot(gs2[4+i])
        ax.text(0,0.8,str(2013+year+1) + ' - ' + str(2013+year+2), fontsize=11, transform=ax.transAxes)
        mixed_similarity = [ 1.-n/N for n in n_random_pairs_mixed_positions[year]]
        plt.hist(mixed_similarity, bins=np.arange(-0.5/N,1+0.5/N,1/N), width=0.8/N, color='lightskyblue')
        adjust_spines(ax,[])
        plt.xlim([-0.5/N,1+0.5/N])

    plt.savefig(os.path.join(work_dir,'mixed_infection_toy_model.png'))
    '''

    '''
    import colorsys
    fig=plt.figure('Allele Frequencies',figsize=(14,8),facecolor='w')
    ax = fig.add_subplot(211)
    def get_colors(ncolors):
        for hue in range(ncolors):
            hue = 1. * hue / ncolors
            col = [int(x) for x in colorsys.hsv_to_rgb(hue, 1.0, 230)]
            yield "#{0:02x}{1:02x}{2:02x}".format(*col)
    color=get_colors(int(N))
    for l in zip(*allele_freqs):
        acolor=next(color)
        plt.plot(l, color=acolor)
    ax.set_ylim([0,1])
    adjust_spines(ax,['left','bottom'])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticks(range(0, n_years+1, 2))
    ax.set_xticklabels(range(2013-n_years+1, 2013+1, 2))
    ax.set_ylabel('Major allele frequency')
    plt.tight_layout()

    Ne=[]
    for y in range(1,len(allele_freqs)):
        F=np.zeros(N)
        for m in range(int(N)):
            # Eq 9 of http://www.genetics.org/content/121/2/379.full.pdf
            # but what about the generation time, which is not the 1-year sampling frequency
            F[m]=(allele_freqs[y-1][m]-allele_freqs[y][m])**2 / ((allele_freqs[y-1][m] + allele_freqs[y][m])/2.0)
        Ne.append(1.0/np.mean(F))
    ax = fig.add_subplot(212)
    ax.plot(Ne)
    adjust_spines(ax,['left'])
    ax.set_ylabel('Variance effective population size (?)')
    plt.tight_layout()
    '''

if __name__ == '__main__':
    #plot_simulation('simulations')
    plot_simulation('simulations/test_likelihoods/sim_8',idx=6) #8:[6,8]
    plt.show()
