from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os

from analyzeCarom import scatter_from_multinest_projection

n_live_list = [20000, 5000, 3000, 2000, 1000]
color_list = ['k', 'r', 'b', 'y', 'c']

physics_dir = os.path.join("/Users", "danielsf", "physics")
fig_dir = os.path.join(physics_dir, "Carom", "figures")
multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")

plt.figsize = (30,30)

header_list = []
label_list= []
for n_live, color in zip(n_live_list, color_list):
    data = None
    for i_fig, dim in enumerate(((6, 9), (0,1))):
        plt.subplot(2,1,i_fig + 1)
        
        file_name = os.path.join(multinest_dir,
                                 "gaussianJellyBean_d12_s99_n%d_t1.00e-03.txt" % n_live)
        
        xx, yy, data = scatter_from_multinest_projection(file_name, 12, dim[0], dim[1],
                                                        data=data)

        hh = plt.scatter(xx, yy, color=color, s=7)
        
        if i_fig != 0:
            scatter_name = file_name.replace('.txt', '_carom.sav')
            with open(scatter_name, 'r') as file_handle:
                lines = file_handle.readlines()
            n_calls = len(lines)-1
            header_list.append(hh)
            label_list.append('n_live %.2e; n_calls %.3e' % (n_live, n_calls))

plt.subplot(2,1,2)
plt.legend(header_list, label_list, fontsize=10, loc=0)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'multinest_nlive_comparison.png'))
