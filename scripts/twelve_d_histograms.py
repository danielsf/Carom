from __future__ import with_statement
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

import scipy.spatial as spatial

from analyzeCarom import  load_dalex_data, load_multinest_data, make_histogram

if __name__ == "__main__":

    seed_list = [626, 694, 762, 1068, 6475]
    limit_list = [400000, 400000, 700000, 400000, 300000]
    title_list = ['(a)', '(b)', '(c)', '(d)', '(e)']
    color_list = ['r', 'b', 'g', 'm', 'k']

    #seed_list = [626, 694]
    #limit_list = [400000, 400000]
    #title_list = ['(a)', '(b)']
    #color_list = ['r', 'b']


    delta_chi = 21.03
    full_dim = 12

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    dalex_dir = os.path.join(physics_dir, "Carom")

    control_dir = os.path.join(dalex_dir, "controls", "draft_160907")
    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    data_dir = os.path.join(dalex_dir, "output", "draft_161022")

    running_min = {}

    plt.figsize = (30, 30)
    for seed, limit, color, title in zip(seed_list, limit_list, color_list, title_list):
        dalex_name = "jellyBean_d12_s%d_output.sav" % seed

        data = load_dalex_data(os.path.join(data_dir, dalex_name), full_dim)
        data = data[:limit]
        r_min = np.minimum.accumulate(data['chisq'])
        running_min[seed] = r_min
        xx, dx = make_histogram(data['chisq'], 1.0, 150.0, cumulative=True)
        plt.plot(xx, dx, color=color, label=title)
    plt.xlabel('$\chi^2$', fontsize=15)
    plt.ylabel('dN/d$\chi^2$', fontsize=15)
    plt.xlim((90.0, 150.0))
    plt.legend(fontsize=20, loc=0)
    plt.savefig(os.path.join(fig_dir, 'dalex_histogram.png'))
    plt.close()
    
    plt.figsize = (30,30)
    for seed, limit, color, title in zip(seed_list, limit_list, color_list, title_list):
        label = title+'; final $\chi^2_{min}=%.2f$'% running_min[seed][-1]
        plt.plot(range(1,len(running_min[seed])+1), running_min[seed], label=label, color=color)
    plt.xlabel('calls to $\chi^2$', fontsize=15)
    plt.ylabel('$\chi^2_{min}$', fontsize=15)
    plt.ylim((90.0, 150.0))
    plt.xlim((1.0e4, 7.0e5))
    plt.legend(fontsize=15, loc=0)
    plt.savefig(os.path.join(fig_dir, 'dalex_minimums.png'))
    plt.close()

    plt.figsize = (30, 30)
    for seed, limit, color, title in zip(seed_list, limit_list, color_list, title_list):
        multinest_name = "gaussianJellyBean_d12_s%d_n300_t1.00e-03.txt" % seed

        data = load_multinest_data(os.path.join(multinest_dir, multinest_name), full_dim)

        xx, dx = make_histogram(data['chisq'], 1.0, 150.0, cumulative=True)
        plt.plot(xx, dx, color=color, label=title)
    plt.xlabel('$\chi^2$', fontsize=15)
    plt.ylabel('dN/d$\chi^2$', fontsize=15)
    plt.xlim((90.0, 150.0))
    plt.legend(fontsize=20, loc=0)
    plt.savefig(os.path.join(fig_dir, 'multinest_histogram.png'))
    plt.close()
