import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from analyzeCarom import make_histogram

if __name__ == "__main__":

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    dalex_dir = os.path.join(physics_dir, "Carom", "output",
                            "workspace")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

    seed_list = [90, 452]
    title_list = ['fast convergence', 'slow convergence']
    cutoff_list = [200000, 400000]

    dtype_list = []
    for ii in range(12):
        dtype_list.append(("p%d" % ii, float))
    dtype_list.append(("chisq", float))
    dtype_list.append(("j1", int))
    dtype_list.append(("j2", int))
    dtype_list.append(("log", int))

    dtype = np.dtype(dtype_list)

    plt.figsize = (30, 30)
    for i_fig, (seed, title, cutoff) in enumerate(zip(seed_list, title_list, cutoff_list)):
        plt.subplot(1, 2, i_fig+1)
        dalex_name = os.path.join(dalex_dir, 'jellyBean_d12_s%d_output.sav' % seed)
        data = np.genfromtxt(dalex_name, dtype=dtype)
        time_list = [50000, 100000] + range(200000, cutoff+50000, 100000)

        header_list = []
        label_list = []
        dchi = 0.5
        chi_cutoff = 130.0
        ymax = 0
        time_list.reverse()
        for time in time_list:
            timed_data = data[:time]
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
        #plt.title(title)
        plt.xlabel('$\chi^2$', fontsize=20)
        plt.ylabel('$dN/d\chi^2$', fontsize=20)
        plt.xlim(data['chisq'].min(),chi_cutoff)
        plt.ylim(0, ymax*1.6)
        print 'ymax ',ymax
        print ''


    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'dalex_histograms.png'))
