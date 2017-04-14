import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from analyzeCarom import make_histogram

if __name__ == "__main__":

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

    n_live = 300

    dtype_list = []
    for ii in range(12):
        dtype_list.append(("p%d" % ii, float))
    dtype_list.append(("chisq", float))
    dtype_list.append(("j1", int))
    dtype_list.append(("j2", int))
    dtype_list.append(("log", int))

    dtype = np.dtype(dtype_list)

    plt.figsize = (30, 30)

    multinest_name = os.path.join(multinest_dir, 'gaussianJellyBean_d12_s99_n%d_t1.00e-03_carom.sav' % n_live)
    data = np.genfromtxt(multinest_name, dtype=dtype)
    cutoff=len(data)
    print('cutoff ',cutoff)
    time_list = [50000, 100000] + range(200000, cutoff, 100000) + [cutoff]

    header_list = []
    label_list = []
    dchi = 0.5
    chi_cutoff = 250.0
    ymax = 0
    time_list.reverse()
    for time in time_list:
        timed_data = data[:time]
        if timed_data['chisq'].min()<chi_cutoff:
            xx, yy = make_histogram(timed_data['chisq'], dchi, chi_cutoff, cumulative=False)
            yy = yy/dchi
            if yy.max()>ymax:
                ymax = yy.max()
            hh, = plt.plot(xx, yy)
            header_list.append(hh)
            label_list.append('%.1e calls; $\chi^2_{min}=%.2f$' % (time, timed_data['chisq'].min()))
            large_dexes = np.where(timed_data['chisq']>chi_cutoff)[0]
            frac = float(len(large_dexes))/float(len(timed_data))
            #print '    frac offending %.6g' % frac
            print '    n offending ',len(large_dexes),' of ',len(timed_data)
            #print '    min ',timed_data['chisq'].min()
            #print ''

    header_list.reverse()
    label_list.reverse()
    plt.legend(header_list, label_list, fontsize=10, loc=0)
    plt.xlabel('$\chi^2$', fontsize=10)
    plt.ylabel('$dN/d\chi^2$', fontsize=10)
    plt.xlim(data['chisq'].min(),chi_cutoff)
    plt.ylim(0, ymax*1.6)
    print 'ymax ',ymax
    print ''
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'multinest_histograms.png'))
