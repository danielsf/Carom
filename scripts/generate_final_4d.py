from __future__ import with_statement
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np
import os

from analyzeCarom import (scatter_from_multinest_projection,
                          scatter_from_carom)

if __name__ == "__main__":

    delta_chi = 9.49

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    dalex_dir = os.path.join(physics_dir, "Carom")

    control_dir = os.path.join(dalex_dir, "controls", "draft_160907")
    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    data_dir = os.path.join(dalex_dir, "output", "workspace")

    n_live = 800

    multinest_file_name = os.path.join(multinest_dir,
                                       "integrableJellyBean_d4_s99_n%d_t1.00e-03.txt" % n_live)

    multinest_scatter_name = os.path.join(multinest_dir,
                                          "integrableJellyBean_d4_s99_n%d_t1.00e-03_carom.sav" % n_live)


    dalex_name = os.path.join(data_dir, "jellyBean_d4_s112_output.sav")

    dim_list = ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3))
    control_dtype = np.dtype([('x', float), ('y', float)])

    plt.figsize = (30, 30)

    m_data = None
    d_data = None

    for i_fig, dim in enumerate(dim_list):
        control_file_name = os.path.join(control_dir,
                                         "gentle_integrable_detailed_x2_0.95_%d_%d_bayesianFullD.txt" %
                                         (dim[0], dim[1]))

        control_data = np.genfromtxt(control_file_name, dtype=control_dtype)

        with open(multinest_scatter_name, 'r') as input_file:
            lines =input_file.readlines()
            n_pts = len(lines) - 1

        m_x, m_y, m_data = scatter_from_multinest_projection(multinest_file_name,
                                                             4, dim[0], dim[1],
                                                             data=m_data,
                                                             downsample=True)

        (d_x, d_y, d_chisq_min, d_target,
         d_data) = scatter_from_carom(dalex_name, 4, dim[0], dim[1], delta_chi=9.49,
                                      data=d_data, limit=n_pts)


        plt.subplot(3,2,i_fig+1)
        hh, = plt.plot(control_data['x'], control_data['y'],
                       color='k', linewidth=1)
        label = 'true 95% projected Bayesian contour'
        if i_fig==0:
            legend_handles = [hh]
            legend_labels = [label]

        hh = plt.scatter(m_x, m_y,color='r', s=5)
        label = '95% projected Bayesian limit (MultiNest)'
        if i_fig==0:
            legend_handles.append(hh)
            legend_labels.append(label)

        hh = plt.scatter(d_x, d_y, color='b', s=10)
        label = '$\chi^2<=\chi^2_{min}+9.49$ contour (Dalex)'
        if i_fig==0:
            legend_handles.append(hh)
            legend_labels.append(label)

        plt.xlabel('$\\theta_%d$' % dim[0], fontsize=15)
        plt.ylabel('$\\theta_%d$' % dim[1], fontsize=15)
        xmin = control_data['x'].min()
        xmax = control_data['x'].max()

        ymin = control_data['y'].min()
        ymax = control_data['y'].max()

        for xc in (d_x.min(), d_x.max(), m_x.min(), m_x.max()):
            if xc < xmin:
                xmin = xc
            if xc > xmax:
                xmax = xc

        for yc in (ymax, ymin, d_y.min(), d_y.max(), m_y.min(), m_y.max()):
            if yc < ymin:
                ymin = yc
            if yc > ymax:
                ymax = yc


        dx = xmax-xmin
        dy = ymax-ymin
        plt.xlim((xmin-0.05*dx, xmax+0.05*dx))
        plt.ylim((ymin-0.05*dy, ymax+0.4*dy))

        if i_fig==0:
            plot_label = '$t=%d$ evaluations of $\chi^2$\n$\chi^2_{min}=%.4e$ (Dalex)' % (n_pts, d_chisq_min)
            plt.text(xmin+0.05*dx, ymax, plot_label, fontsize=10)

    plt.tight_layout()
    file_name = os.path.join(fig_dir, "final_4d.eps")
    plt.savefig(file_name)
