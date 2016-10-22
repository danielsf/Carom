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
    data_dir = os.path.join(dalex_dir, "output", "draft_161011")

    multinest_file_list = []
    multinest_scatter_file_list = []
    for ix in (100, 300, 400, 500):
        multinest_file_name = os.path.join(multinest_dir,
                                           "integrableJellyBean_d4_s99_n%d_t1.00e-03.txt" % ix)

        multinest_file_list.append(multinest_file_name)

        multinest_file_name = os.path.join(multinest_dir,
                                           "integrableJellyBean_d4_s99_n%d_t1.00e-03_carom.sav" % ix)

        multinest_scatter_file_list.append(multinest_file_name)

    dim_list = ((0,1), (0,2), (0,3), (1,2), (1,3), (2,3))
    
    control_dtype = np.dtype([('x', float), ('y', float)])

    control_data = {}
    for dim in dim_list:
        control_file_name = os.path.join(control_dir,
                                         "gentle_integrable_detailed_x2_0.95_%d_%d_bayesianFullD.txt" %
                                         (dim[0], dim[1]))

        control_data[dim] = np.genfromtxt(control_file_name, dtype=control_dtype)

    dalex_name = os.path.join(data_dir, "jellyBean_d4_output.sav")
    
    (x_grid, y_grid, dalex_chisq_min, dalex_target,
     dalex_data) = scatter_from_carom(dalex_name, 4, 0, 1, delta_chi=9.49)

    legend_handles = {}
    legend_labels = {}
    for dim in dim_list:
        legend_handles[dim] = []
        legend_labels[dim] = []

    xmin_dict = {}
    xmax_dict = {}
    ymin_dict = {}
    ymax_dict = {}
    plot_labels_dict = {}
    for dim in dim_list:
        plot_labels_dict[dim] = []

    for i_time, (multinest_file, multinest_scatter) in \
    enumerate(zip(multinest_file_list, multinest_scatter_file_list)):

        with open(multinest_scatter, 'r') as input_file:
            lines =input_file.readlines()
            n_pts = len(lines) - 1

        for i_fig, dim in enumerate(dim_list):
            plt.figure(i_fig+1)
            plt.subplot(3,2,i_time+1)
            m_x, m_y, m_data = scatter_from_multinest_projection(multinest_file,
                                      4, dim[0], dim[1])
            
            (d_x, d_y, d_chisq_min, d_target,
             d_data) = scatter_from_carom("junk.txt", 4, dim[0], dim[1], delta_chi=9.49,
                                          data=dalex_data, limit=n_pts)

            hh, = plt.plot(control_data[dim]['x'], control_data[dim]['y'],
                           color='k', linewidth=1)
            label = 'true 95% projected Bayesian contour'
            if label not in legend_labels[dim]:
                legend_handles[dim].append(hh)
                legend_labels[dim].append(label)
            
            hh = plt.scatter(m_x, m_y,color='r', s=5)
            label = '95% projected Bayesian limit (MultiNest)'
            if label not in legend_labels[dim]:
                legend_handles[dim].append(hh)
                legend_labels[dim].append(label)
            
            hh = plt.scatter(d_x, d_y, color='b', s=10)
            label = '$\chi^2<=\chi^2_{min}+9.49$ contour (Dalex)'
            if label not in legend_labels[dim]:
                legend_handles[dim].append(hh)
                legend_labels[dim].append(label)
            
            plt.xlabel('$\\theta_%d$' % dim[0], fontsize=10)
            plt.ylabel('$\\theta_%d$' % dim[1], fontsize=10)
            xmin = control_data[dim]['x'].min()
            xmax = control_data[dim]['x'].max()
            
            ymin = control_data[dim]['y'].min()
            ymax = control_data[dim]['y'].max()
            
            if dim not in xmax_dict:
                xmax_dict[dim] = xmax
                xmin_dict[dim] = xmin
                ymax_dict[dim] = ymax
                ymin_dict[dim] = ymin
            
            for xc in (xmax, xmin, d_x.min(), d_x.max(), m_x.min(), m_x.max()):
                if xc < xmin_dict[dim]:
                    xmin_dict[dim] = xc
                if xc > xmax_dict[dim]:
                    xmax_dict[dim] = xc
            
            for yc in (ymax, ymin, d_y.min(), d_y.max(), m_y.min(), m_y.max()):
                if yc < ymin_dict[dim]:
                    ymin_dict[dim] = yc
                if yc > ymax_dict[dim]:
                    ymax_dict[dim] = yc
            
            
            label = '$t=%d$ evaluations of $\chi^2$\n$\chi^2_{min}=%.4e$ (Dalex)' % (n_pts, d_chisq_min)
            plot_labels_dict[dim].append(label)
            
                

    for i_fig, dim in enumerate(dim_list):
        plt.figure(i_fig+1)
        for i_time in range(1,len(multinest_file_list)+1):
            plt.subplot(3,2,i_time)
            dx = xmax_dict[dim] - xmin_dict[dim]
            plt.xlim((xmin_dict[dim]-0.05*dx, xmax_dict[dim]+0.05*dx))
            dy = ymax_dict[dim] - ymin_dict[dim]
            plt.ylim((ymin_dict[dim]-0.05*dy, ymax_dict[dim]+0.4*dy))
            
            plt.text(xmin_dict[dim]+0.05*dx, ymax_dict[dim],
                     plot_labels_dict[dim][i_time-1], fontsize=10)
        
        #plt.subplot(3,2,i_time+1)    
        plt.legend(legend_handles[dim], legend_labels[dim], fontsize=10,
                   bbox_to_anchor=(-0.05,0.0), loc=2)
        plt.tight_layout()
        file_name = os.path.join(fig_dir, "time_lapse_%d_%d.png" % (dim[0], dim[1]))
        plt.savefig(file_name)
