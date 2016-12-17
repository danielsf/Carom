from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os

from analyzeCarom import scatter_from_multinest_projection
from analyzeCarom import scatter_from_carom


physics_dir = os.path.join("/Users", "danielsf", "physics")
carom_dir = os.path.join(physics_dir, "Carom")
fig_dir = os.path.join(carom_dir, "figures")
dalex_dir = os.path.join(carom_dir, "output", "draft_161117")
multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")

dalex_limit = 250000
dalex_name = os.path.join(dalex_dir, "lump_d12_s13_output.sav")

n_list = [10000, 5000, 3000, 2000, 1000]
color_list = ['k', 'r', 'b', 'y', 'c']

plt.figsize = (30, 30)

m_header_list = []
m_label_list = []

xmin = {}
xmax = {}
ymin = {}
ymax = {}

color_map = plt.cm.gist_ncar

for i_set, n_live in enumerate(n_list):
    multinest_name = os.path.join(multinest_dir,
                                  "nonGaussianLump_d12_s99_n%d_t1.00e-03.txt"
                                  % (n_live))

    data = None
    for i_fig, dim in enumerate(((6, 9), (0, 3))):
        plt.subplot(2,2,i_fig*2+1)
        xx, yy, data = scatter_from_multinest_projection(multinest_name, 12, dim[0], dim[1],
                                                         data=data)

        if dim not in xmin:
            xmin[dim] = xx.min()
            xmax[dim] = xx.max()
            ymin[dim] = yy.min()
            ymax[dim] = yy.max()
        else:
            if xx.min()<xmin[dim]:
                xmin[dim] = xx.min()
            if xx.max()>xmax[dim]:
                xmax[dim] = xx.max()
            if yy.min()<ymin[dim]:
                ymin[dim]=yy.min()
            if yy.max()>ymax[dim]:
                ymax[dim]=yy.max()

        color= color_map(i_set*50)
        hh = plt.scatter(xx, yy, color=color, s=7)
        if i_fig != 0:
            scatter_name = multinest_name.replace('.txt', '_carom.sav')
            with open(scatter_name, 'r') as file_handle:
                lines = file_handle.readlines()
            n_calls = len(lines)-1
            m_header_list.append(hh)
            m_label_list.append('%.1e live; %.2e calls' % (n_live, n_calls))


multinest_name = os.path.join(multinest_dir,
                                  "nonGaussianLump_d12_s99_n10000_t1.00e-03.txt")

control_data = None
dalex_data = None

for i_fig, dim in enumerate(((6, 9), (0, 3))):
    plt.subplot(2,2,i_fig*2+2)
    d_x, d_y, d_min, d_targ, dalex_data = scatter_from_carom(dalex_name, 12, dim[0], dim[1],
                                                             data=dalex_data, delta_chi=21.03,
                                                             limit=dalex_limit)

    m_x, m_y, control_data = scatter_from_multinest_projection(multinest_name, 12, dim[0], dim[1],
                                                               data=control_data)

    plt.scatter(m_x, m_y, color='k', s=7)
    plt.scatter(d_x, d_y, color='r', s=7, marker='+')

    for xx in (d_x.min(), d_x.max(), m_x.min(), m_x.max()):
        if xx<xmin[dim]:
            xmin[dim]=xx
        if xx>xmax[dim]:
            xmax[dim]=xx

    for yy in (d_y.min(), d_y.max(), m_y.min(), m_y.max()):
        if yy<ymin[dim]:
            ymin[dim]=yy
        if yy>ymax[dim]:
            ymax[dim]=yy


for i_fig in (1,2):
    plt.subplot(2,2,i_fig)
    dim=(6,9)
    dx=xmax[dim]-xmin[dim]
    dy=ymax[dim]-ymin[dim]
    plt.xlim((xmin[dim]-0.05*dx,xmax[dim]+0.05*dx))
    plt.ylim((ymin[dim]-0.05*dx,ymax[dim]+0.05*dy))

for i_fig in (3,4):
    plt.subplot(2,2,i_fig)
    dim=(0,3)
    dx=xmax[dim]-xmin[dim]
    dy=ymax[dim]-ymin[dim]
    plt.xlim((xmin[dim]-0.05*dx,xmax[dim]+0.1*dx))
    plt.ylim((ymin[dim]-0.05*dy,ymax[dim]+0.1*dy))
    if i_fig==3:
        plt.legend(m_header_list, m_label_list, fontsize=7, loc=0)

plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "lump_n_live_comparison.png"))
