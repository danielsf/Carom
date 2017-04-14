from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os

from analyzeCarom import scatter_from_multinest_projection

n_live_list = [50000, 20000, 5000, 3000, 2000]
color_list = ['k', 'r', 'b', 'y', 'c']

physics_dir = os.path.join("/Users", "danielsf", "physics")
fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")

plt.figsize = (30,30)

header_list = []
label_list= []

xmin = 3.0e30
xmax = -3.0e30
ymin = 3.0e30
ymax = -3.0e30

for i_set, n_live in enumerate(n_live_list):
    print 'plotting ',i_set
    data = None
    for i_fig, dim in enumerate(((6, 9), (0,1))):
        plt.subplot(2,1,i_fig + 1)

        file_name = os.path.join(multinest_dir,
                                 "gaussianJellyBean_d12_s99_n%d_t1.00e-03.txt" % n_live)

        xx, yy, data = scatter_from_multinest_projection(file_name, 12, dim[0], dim[1],
                                                        data=data,
                                                        downsample=0.02)

        color = plt.cm.gist_ncar(i_set*50)
        hh = plt.scatter(xx, yy, color=[color]*len(xx), s=7)

        plt.xlabel('$\\theta_%d$' % dim[0], fontsize=15)
        plt.ylabel('$\\theta_%d$' % dim[1], fontsize=15)

        if i_fig != 0:
            if xx.min()<xmin:
                xmin = xx.min()
            if xx.max()>xmax:
                xmax = xx.max()
            if yy.min()<ymin:
                ymin = yy.min()
            if yy.max()>ymax:
                ymax = yy.max()
            scatter_name = file_name.replace('.txt', '_carom.sav')
            with open(scatter_name, 'r') as file_handle:
                lines = file_handle.readlines()
            n_calls = len(lines)-1
            header_list.append(hh)
            label_list.append('n_live %.2e; n_calls %.3e' % (n_live, n_calls))

plt.subplot(2,1,2)

dx = xmax-xmin
dy = ymax-ymin

plt.xlim((xmin, xmax+0.4*dx))
plt.ylim((ymin, ymax+0.4*dx))
plt.legend(header_list, label_list, fontsize=10, loc=0)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'jellyBean_nlive_comparison.eps'))
