import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os

from analyzeCarom import scatter_from_multinest_projection
from analyzeCarom import scatter_from_carom

import sys

if __name__ == "__main__":

    dim = (int(sys.argv[1]), int(sys.argv[2]))

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    dalex_dir = os.path.join(physics_dir, "Carom")

    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    data_dir = os.path.join(dalex_dir, "output", "draft_161111")

    dalex_name = os.path.join(data_dir, "lump_output.sav")
    multinest_name = os.path.join(multinest_dir,
                                 "nonGaussianLump_d12_s99_n300_t1.00e-03.txt")

    multinest_carom_name = os.path.join(multinest_dir,
                                        "nonGaussianLump_d12_s99_n300_t1.00e-03_carom.sav")

    plt.figsize = (30,30)
    time_list = [75000, 100000, 125000, 150000, 200000]
    delta_chisq = 21.03
    full_dim = 12

    d_data = None
    m_data = None
    xmax = -1.0e30
    xmin = 1.0e30
    ymax = -1.0e30
    ymin = 1.0e30
    mult_min_list = []
    dalex_min_list = []
    m_h = None
    for i_fig, limit in enumerate(time_list):
        plt.subplot(3, 2, i_fig+1)

        if i_fig == len(time_list)-1:
            m_x_true, m_y_true, _data = scatter_from_multinest_projection(multinest_name,
                                                                full_dim, dim[0], dim[1])

            for xx in (m_x_true.min(), m_x_true.max()):
                if xx>xmax:
                    xmax=xx
                if xx<xmin:
                    xmin=xx
            for yy in (m_y_true.min(), m_y_true.max()):
                if yy>ymax:
                    ymax=yy
                if yy<ymin:
                    ymin=yy


        (d_x, d_y, chisq_min_dalex, target_dalex,
         d_data) = scatter_from_carom(dalex_name, full_dim, dim[0], dim[1],
                                      data=d_data,
                                      delta_chi=delta_chisq, limit=limit)

        (m_x, m_y, chisq_min_mult, target_mult,
         m_data) = scatter_from_carom(multinest_carom_name, full_dim, dim[0], dim[1],
                                      data=m_data,
                                      delta_chi=delta_chisq, limit=limit)


        time_list[len(time_list)-1] = len(m_data)
        mult_min_list.append(chisq_min_mult)
        dalex_min_list.append(chisq_min_dalex)

        for xx in (m_x.min(), m_x.max(), d_x.min(), d_x.max()):
            if xx>xmax:
                xmax=xx
            if xx<xmin:
                xmin=xx

        for yy in (m_y.min(), m_y.max(), d_y.min(), d_y.max()):
            if yy>ymax:
                ymax=yy
            if yy<ymin:
                ymin=yy

        d_h = plt.scatter(d_x, d_y, color='r', s=7)

        m_color = 'b'
        if i_fig == len(time_list)-1:
             m_true_h = plt.scatter(m_x_true, m_y_true, color='b', s=7)
             m_color = 'c'

        _m_h = plt.scatter(m_x, m_y, color=m_color, marker='+', s=20)
        if m_h is None:
            m_h = _m_h


    for i_fig in range(len(time_list)):
        plt.subplot(3,2,i_fig+1)
        dx = xmax-xmin
        dy = ymax-ymin
        plt.xlim((xmin-0.05*dx, xmax+0.05*dx))
        plt.ylim((ymin-0.05*dy, ymax+0.5*dy))
        if i_fig==0:
            plt.xlabel('$\\theta_%d$' % dim[0], fontsize=15)
            plt.ylabel('$\\theta_%d$' % dim[1], fontsize=15)
        msg = ('MultiNest $\chi^2_{min}=%.2f$\nDale$\chi$ $\chi^2_{min}=%.2f$\n$\chi^2$ evaluations: %.2e'
               % (mult_min_list[i_fig], dalex_min_list[i_fig], time_list[i_fig]))

        if dim[0]==0 and dim[1]==3:
            xorig = xmin+0.4*dx
            yorig = ymax-0.25*dy
        else:
            xorig = xmin+0.01*dx
            yorig = ymax-0.25*dy

        plt.text(xorig, yorig, msg, fontsize=10)

    plt.legend([d_h, m_h, m_true_h],
               ['Dale$\chi$ ($\chi^2<=\chi^2_{min}+%.2f$)' % delta_chisq,
                'MultiNest ($\chi^2<=\chi^2_{min}+%.2f$)' % delta_chisq,
                'MultiNest (Bayesian, projected)'],
                fontsize=10,
                bbox_to_anchor=(1.05,1),
                loc=2)

    fig_name = 'lump_%d_%d.png' % (dim[0], dim[1])
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, fig_name))
    plt.close()
    print xmin,xmax,ymin,ymax
