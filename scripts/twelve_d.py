from __future__ import with_statement
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

from analyzeCarom import scatter_from_multinest_projection, scatter_from_carom
from analyzeCarom import scatter_from_multinest_marginalized

if __name__ == "__main__":

    seed_list = [932, 99, 66, 125, 6475]
    title_list = ['(a)', '(b)', '(c)', '(d)', '(e)']
    limit_list = [500000, 400000, 300000, 200000]
    color_list = ['k', 'y', 'g', 'b']

    delta_chi = 21.03
    full_dim = 12

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    dalex_dir = os.path.join(physics_dir, "Carom")

    control_dir = os.path.join(dalex_dir, "controls", "draft_160907")
    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    data_dir = os.path.join(dalex_dir, "output", "draft_160907")

    dim_list = [(6,9), (0,1)]

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

        for i_seed, (title, seed) in enumerate(zip(title_list, seed_list)):
            dalex_name = "jellyBean_d12_s%d_output.sav" % seed
            multinest_name = "gaussianJellyBean_d12_s%d_n300_t1.00e-03.txt" % seed

            scatter_name = "gaussianJellyBean_d12_s%d_n300_t1.00e-03_carom.sav" % seed
            with open(os.path.join(multinest_dir, scatter_name), 'r') as input_file:
                lines = input_file.readlines()
                n_mult = len(lines)


            m_x, m_y, m_data = scatter_from_multinest_projection(
                                 os.path.join(multinest_dir, multinest_name),
                                 full_dim, dim[0], dim[1],
                                 data=m_data_dict[seed])

            m_data_dict[seed] = m_data

            plt.subplot(3,2,i_seed+1)
            plt.title(title, fontsize=20)
            m_h = plt.scatter(m_x, m_y, color='r', s=7)

            d_h_list = []
            d_label_list = []

            for limit, color in zip(limit_list, color_list):

                (d_x, d_y, chisq_min, target,
                 d_data) = scatter_from_carom(os.path.join(data_dir, dalex_name),
                                          full_dim, dim[0], dim[1], delta_chi=delta_chi,
                                          data=d_data_dict[seed],
                                          limit=limit)

                d_data_dict[seed] = d_data

                d_h = plt.scatter(d_x, d_y, color=color, s=3)
                d_h_list.append(d_h)
                d_label_list.append('Dale$\chi$: %.2e calls' % limit)

                if xmax is None or d_x.max()>xmax:
                    xmax=d_x.max()
                if xmin is None or d_x.min()<xmin:
                    xmin=d_x.min()
                if ymax is None or d_y.max()>ymax:
                    ymax=d_y.max()
                if ymin is None or d_y.min()<ymin:
                    ymin=d_y.min()

            text = ('MultiNest: %.2e $\chi^2$ calls\n' % (n_mult))
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

        for i_seed in range(5):
            plt.subplot(3,2,i_seed+1)
            dx=xmax-xmin
            dy=ymax-ymin
            plt.xlim((xmin-0.05*dx, xmax+0.05*dx))
            plt.ylim((ymin-0.05*dy, ymax+0.6*dy))
            plt.xlabel('$\\theta_%d$' % dim[0], fontsize=20)
            plt.ylabel('$\\theta_%d$' % dim[1], fontsize=20)
            plt.text(xmin, ymax+0.01*dy, text_list[i_seed], fontsize=10)

        plt.legend([m_h] + d_h_list,
                   ['MultiNest'] + d_label_list,
                   fontsize=10,
                   bbox_to_anchor=(1.05,1),
                   loc=2)

        fig_name = 'compare_d12_%d_%d.png' % (dim[0], dim[1])
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, fig_name))
        plt.close()
