from __future__ import with_statement
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

import scipy.spatial as spatial

from analyzeCarom import scatter_from_multinest_projection, scatter_from_carom
from analyzeCarom import scatter_from_multinest_marginalized


if __name__ == "__main__":

    seed_list = [83, 13, 48, 33, 62]
    limit_list = [250000, 150000, 250000, 250000, 250000]
    title_list = ['(a)', '(b)', '(c)', '(d)', '(e)']
    color_list = ['r', 'y', 'g', 'c']

    delta_chi = 21.03
    full_dim = 12

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    dalex_dir = os.path.join(physics_dir, "Carom")

    control_dir = os.path.join(dalex_dir, "controls", "draft_160907")
    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    data_dir = os.path.join(dalex_dir, "output", "draft_161117")

    dim_list = [(6,9), (0,3)]

    d_data_dict = {}
    m_data_dict = {}
    for seed in seed_list:
        d_data_dict[seed] = None
        m_data_dict[seed] = None


    for dim in dim_list:
        plt.figsize = (30,30)

        text_list = []
        xmax = None
        xmin = None
        ymax = None
        ymin = None

        for i_seed, (title, seed, limit) in enumerate(zip(title_list, seed_list, limit_list)):
            t_start = time.time()
            dalex_name = "lump_d12_s%d_output.sav" % seed
            multinest_name = "nonGaussianLump_d12_s%d_n300_t1.00e-03.txt" % seed

            scatter_name = "nonGaussianLump_d12_s%d_n300_t1.00e-03_carom.sav" % seed
            with open(os.path.join(multinest_dir, scatter_name), 'r') as input_file:
                lines = input_file.readlines()
                n_mult = len(lines)

            plt.subplot(3,2,i_seed+1)
            plt.title(title, fontsize=7)

            d_h_list = []
            d_label_list = []

            (d_x, d_y, chisq_min, target,
             d_data) = scatter_from_carom(os.path.join(data_dir, dalex_name),
                                          full_dim, dim[0], dim[1], delta_chi=delta_chi,
                                          data=d_data_dict[seed],
                                          limit=limit)

            d_data_dict[seed] = d_data

            if xmax is None or d_x.max()>xmax:
                xmax=d_x.max()
            if xmin is None or d_x.min()<xmin:
                xmin=d_x.min()
            if ymax is None or d_y.max()>ymax:
                ymax=d_y.max()
            if ymin is None or d_y.min()<ymin:
                ymin=d_y.min()

            m_x, m_y, m_data = scatter_from_multinest_projection(
                                 os.path.join(multinest_dir, multinest_name),
                                 full_dim, dim[0], dim[1],
                                 data=m_data_dict[seed])

            m_data_dict[seed] = m_data

            m_h = plt.scatter(m_x, m_y, color='b', s=7)


            d_h = plt.scatter(d_x, d_y, color='r', s=7, marker='+')
            d_h_list.append(d_h)
            d_label_list.append('Dale$\chi$; $\chi^2<=\chi^2_{min}+21.03$')

            text = ('MultiNest: %.2e $\chi^2$ calls\n' % n_mult
                    + 'Dale$\chi$: %.2e $\chi^2$ calls; $\chi^2_{min}=%.2f$' % (limit, chisq_min))
            text_list.append(text)

            for xx in (m_x.min(), m_x.max(), d_x.min(), d_x.max()):
                if xmax is None or xx>xmax:
                    xmax=xx
                if xmin is None or xx<xmin:
                    xmin=xx

            for yy in (m_y.min(), m_y.max(), d_y.min(), d_y.max()):
                if ymax is None or yy>ymax:
                    ymax=yy
                if ymin is None or yy<ymin:
                    ymin=yy
            print "one seed took ",time.time()-t_start

        for i_seed in range(len(seed_list)):
            plt.subplot(3,2,i_seed+1)
            dx=xmax-xmin
            dy=ymax-ymin
            plt.xlim((xmin-0.05*dx, xmax+0.05*dx))
            plt.ylim((ymin-0.05*dy, ymax+dy))
            if i_seed==0:
                plt.xlabel('$\\theta_%d$' % dim[0], fontsize=15)
                plt.ylabel('$\\theta_%d$' % dim[1], fontsize=15)
            plt.text(xmin, ymax+0.1*dy, text_list[i_seed], fontsize=10)

        plt.legend([m_h] + d_h_list,
                   ['MultiNest'] + d_label_list,
                   fontsize=10,
                   bbox_to_anchor=(1.05,1),
                   loc=2)

        fig_name = 'compare_d12_lump_%d_%d.png' % (dim[0], dim[1])
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, fig_name))
        plt.close()
