from __future__ import with_statement
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from analyzeDalex import scatter_from_multinest_projection, scatter_from_dalex
from analyzeDalex import scatter_from_multinest_marginalized


if __name__ == "__main__":

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    dalex_dir = os.path.join(physics_dir, "Carom", "output", "workspace")

    delta_chi = 21.03
    full_dim = 12
    nlive = 20000
    n_dud = 300
    seed = 90

    multinest_dud_file = os.path.join(multinest_dir,
                                      "gaussianJellyBean_d12_s99_n%d_t1.00e-03.txt" % n_dud)

    multinest_dud_scatter = os.path.join(multinest_dir,
                                         "gaussianJellyBean_d12_s99_n%d_t1.00e-03_carom.sav" % n_dud)


    multinest_file = os.path.join(multinest_dir,
                                  "gaussianJellyBean_d12_s99_n%d_t1.00e-03.txt" % nlive)

    multinest_scatter = os.path.join(multinest_dir,
                                     "gaussianJellyBean_d12_s99_n%d_t1.00e-03_carom.sav" % nlive)


    with open(multinest_scatter, 'r') as file_handle:
        scatter_lines = file_handle.readlines()
    n_mult_calls = len(scatter_lines)-1
    del scatter_lines

    with open(multinest_dud_scatter, 'r') as file_handle:
        scatter_lines = file_handle.readlines()
    n_mult_dud_calls = len(scatter_lines)-1
    del scatter_lines

    dalex_file = os.path.join(dalex_dir,'jellyBean_d12_s%d_output.sav' % seed)

    m_data = None
    md_data = None
    d_data = None

    time_list = [150000, 200000, 300000, 400000]

    for ix, iy in zip((0, 6), (1, 9)):

        m_x, m_y, m_data = scatter_from_multinest_projection(multinest_file,
                                                             full_dim, ix, iy, data=m_data,
                                                             downsample=0.02)

        md_x, md_y, md_data = scatter_from_multinest_projection(multinest_dud_file,
                                                                full_dim, ix, iy, data=md_data,
                                                                downsample=0.02)


        print('len %d %d' % (len(m_x), len(md_x)))

        xmin = m_x.min()
        xmax = m_x.max()
        ymin = m_y.min()
        ymax = m_y.max()

        plt.figsize = (30, 30)
        header_list = []
        label_list = []
        min_dict = {}
        for i_fig, limit in enumerate(time_list):

            plt.subplot(2,2,i_fig+1)
            m_hh = plt.scatter(m_x, m_y, color='k', s=5)
            m_d_hh = plt.scatter(md_x, md_y, color='c', s=5)
            if i_fig == len(time_list)-1:
                header_list.append(m_hh)
                label_list.append('MultiNest; %.2e $\chi^2$ calls' % n_mult_calls)
                header_list.append(m_d_hh)
                label_list.append('MultiNest; %.2e $\chi^2$ calls' % n_mult_dud_calls)

            (d_x, d_y, d_min,
             d_target, d_data) = scatter_from_dalex(dalex_file, full_dim, ix, iy,
                                                    delta_chi=delta_chi, data=d_data,
                                                    limit=limit)

            min_dict[limit] = d_min

            if d_x.min()<xmin:
                xmin=d_x.min()
            if d_x.max()>xmax:
                xmax=d_x.max()
            if d_y.min()<ymin:
                ymin=d_y.min()
            if d_y.max()>ymax:
                ymax=d_y.max()

            d_hh = plt.scatter(d_x, d_y, color='r', marker='x', s=7)
            if i_fig == len(time_list)-1:
                header_list.append(d_hh)
                label_list.append('Dale$\chi$')

        dx=xmax-xmin
        dy=ymax-ymin

        xmin -= 0.2*dx
        xmax += 0.1*dx
        ymin -= 0.4*dy
        ymax += 0.35*dy

        for i_fig, limit in enumerate(time_list):
            plt.subplot(2,2,i_fig+1)

            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            plt.text(xmin+0.05*dx, ymax-0.3*dy,
                     'Dale$\chi$ calls to $\chi^2$: %.2e\n' % limit
                     + 'Dale$\chi$ $\chi^2$ min: %.3e' % min_dict[limit],
                     fontsize=10)

            if i_fig==len(time_list)-1:
                header_list.reverse()
                label_list.reverse()
                plt.legend(header_list, label_list, fontsize=10, loc=3,
                           frameon=False)
            elif i_fig==0:
                plt.xlabel('$\\theta_%d$' % ix, fontsize=15)
                plt.ylabel('$\\theta_%d$' % iy, fontsize=15)

        plt.tight_layout()

        if seed==90 and ix==0 and iy==1:
            file_name = os.path.join(fig_dir, 'figure_11.png')
        elif seed==90 and ix==6 and iy==9:
            file_name = os.path.join(fig_dir, 'figure_10.png')
        else:
            file_name = os.path.join(fig_dir,
                                     'jellyBean_time_series_s%d_%d_%d.eps' % (seed, ix, iy))

        plt.savefig(file_name)
        plt.close()

