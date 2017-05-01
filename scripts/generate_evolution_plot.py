import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyzeDalex import scatter_from_multinest_projection
from analyzeDalex import scatter_from_dalex

import os

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--x', type=int, default=0)
parser.add_argument('--y', type=int, default=1)
args = parser.parse_args()

seed_list = [452, 19, 642, 5345, 7742, 831]
limit_list = [800000, 600000, 400000]
color_list = ['c', 'y', 'r']

physics_dir = os.path.join('/Users', 'danielsf', 'physics')

out_dir = os.path.join(physics_dir, 'Carom_drafts', 'figures')

dalex_dir = os.path.join(physics_dir, 'Carom', 'output', 'followup')

control_file = os.path.join(physics_dir, 'MultiNest_v3.9',
                            'chains',
                            'gaussianJellyBean_d12_s99_n50000_t1.00e-03.txt')

(m_x,
 m_y, 
 m_data) = scatter_from_multinest_projection(control_file, 12, args.x, args.y,
                                             downsample=0.02)

plt.figsize = (30,30)

for i_fig, seed in enumerate(seed_list):
    dalex_file = os.path.join(dalex_dir, 'jellyBean_d12_s%d_output.sav' % seed)
    plt.subplot(3,2,i_fig+1)
    
    header_list = []
    label_list = []
    
    hh_c = plt.scatter(m_x, m_y, color='k', marker='o', s=5)
    if i_fig == 0:
        plt.legend([hh_c], ['MultiNest after %.2e samples' % 35000000],
                   fontsize=10, loc=0)
    
    data = None
    for i_c, (lim, color) in enumerate(zip(limit_list, color_list)):
        (d_x, d_y,
         d_min, d_target,
         data) = scatter_from_dalex(dalex_file, 12, args.x, args.y,
                                    delta_chi=21.03, data=data, limit=lim)

        hh = plt.scatter(d_x, d_y, color=color, s=5)
        if i_fig == 3-i_c:
            plt.legend([hh], ['Dale$\chi$ after %.2e samples' % lim],
                       fontsize=10, loc=0)

    if i_fig==0:
        plt.xlabel('$\\theta_%d$' % args.x)
        plt.ylabel('$\\theta_%d$' % args.y)

    plt.xlim(-20, 15)
    plt.ylim(0,90)

plt.tight_layout()
out_name = os.path.join(out_dir,'jellyBean_evolution.eps')
plt.savefig(out_name)
