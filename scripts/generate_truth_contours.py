import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import os
import numpy as np


def order_x_y(x_in, y_in):
    available = range(len(x_in))
    available.pop(0)

    x_out = np.zeros(len(x_in))
    y_out = np.zeros(len(y_in))

    x_out[0] = x_in[0]
    y_out[0] = y_in[0]

    steps = np.zeros(len(x_in))
    n_chosen = 1

    while n_chosen != len(x_in):
        distance = [(x_out[n_chosen-1]-x_in[ii])**2 +
                    (y_out[n_chosen-1]-y_in[ii])**2
                    for ii in available]

        distance = np.array(distance)
        chosen_dex = np.argmin(distance)
        chosen_distance = distance[chosen_dex]

        if n_chosen>2:
            mean_step = np.mean(steps[:len(x_out)])
            std_step = np.std(steps[:len(x_out)])

        if n_chosen<=2 or chosen_distance<1.0:
            steps[n_chosen] = chosen_distance
            x_out[n_chosen] = x_in[available[chosen_dex]]
            y_out[n_chosen] = y_in[available[chosen_dex]]

        else:
            distance = (np.power(x_out[:n_chosen]-x_in[available[chosen_dex]],2) +
                        np.power(y_out[:n_chosen]-y_in[available[chosen_dex]],2))

            insert_dex = np.argmin(distance)
            for ii in range(n_chosen, insert_dex+1, -1):
                x_out[ii] = x_out[ii-1]
                y_out[ii] = y_out[ii-1]

            x_out[insert_dex+1] = x_in[available[chosen_dex]]
            y_out[insert_dex+1] = y_in[available[chosen_dex]]

            steps[0] = 0.0
            for ii in range(1, n_chosen):
                steps[ii] = (np.power(x_out[ii]-x_out[ii-1],2) +
                             np.power(y_out[ii]-y_out[ii-1],2))

        n_chosen += 1
        available.pop(chosen_dex)

    return x_out, y_out


if __name__ == "__main__":

    control_dir = os.path.join("/Users", "danielsf", "physics",
                               "Carom", "controls", "draft_160907")

    output_dir = os.path.join("/Users", "danielsf", "physics",
                              "Carom_drafts", "figures")

    control_root = "gentle_integrable_detailed_x2"

    dtype = np.dtype([('x', float), ('y', float)])

    plt.figsize = (30,30)

    i_fig = 0

    legend_dict = {}
    legend_dict['frequentistFullDrelative'] = {}
    legend_dict['frequentistFullDrelative'][0.95] = \
    '$\chi^2 <= \chi^2_{min} + 9.49 \, (4.70)$'
    legend_dict['frequentistFullDrelative'][0.68] = \
    '$\chi^2 <= \chi^2_{min} + 4.70$'
    
    legend_dict['frequentistFullD'] = {}
    legend_dict['frequentistFullD'][0.95] = \
    '$\chi^2 <= 124.35$'
    legend_dict['frequentistFullD'][0.68] = \
    '$\chi^2 <= 106.07$'
    
    legend_dict['frequentist2Deff'] = {}
    legend_dict['frequentist2Deff'][0.95] = \
    '$\chi^2_{eff} <= \chi^2_{eff, min} + 6.0 \, (2.28)$' # this is not quite correct
    legend_dict['frequentist2Deff'][0.68] = \
    '$\chi^2_{eff} <= \chi^2_{eff, min} + 2.28$' # this is not quite correct
    
    legend_dict['frequentist2D'] = {}
    legend_dict['frequentist2D'][0.95] = \
    '$\chi^2 <= \chi^2_{min} + 6.0 \, (2.28)$' # this is not quite correct
    legend_dict['frequentist2D'][0.68] = \
    '$\chi^2 <= \chi^2_{min} + 2.28$' # this is not quite correct
    
    
    legend_dict['bayesianFullD'] = {}
    legend_dict['bayesianFullD'][0.95] = \
    '95% (68%) of Bayesian posterior (proj.)'
    legend_dict['bayesianFullD'][0.68] = \
    '68% of Bayesian posterior (proj.)'
    
    legend_dict['bayesian2D'] = {}
    legend_dict['bayesian2D'][0.95] = \
    '95% (68%) of Bayesian posterior (marg.)'
    legend_dict['bayesian2D'][0.68] = \
    '68% of Bayesian posterior (marg.)'


    n_rows = 3
    n_cols = 2
    iplot = 0

    for ix in range(0, 4):
        for iy in range(ix+1,4):
            legend_handles = []
            legend_labels = []
            if ix<iy and iy!=1:
                i_fig += 1
                plt.subplot(n_rows, n_cols, i_fig)
                
                for suffix, color in zip(('bayesian2D', 'frequentist2D',
                                          'bayesianFullD', 'frequentistFullDrelative'
                                          ),
                                          ('r', 'g', 'b', 'm')):
                
                    for pct in (0.68, 0.95):
                    
                        file_name = '%s_%.2f_%d_%d_%s.txt' % \
                        (control_root, pct, ix, iy, suffix)
                        
                        data = np.genfromtxt(os.path.join(control_dir, file_name),
                                             dtype=dtype)
                
                        label = legend_dict[suffix][pct]
                        if pct>0.7:
                            linewidth=2
                        else:
                            linewidth=0.5
                        
                        #print suffix,linestyle,file_name,len(data)
                        hh, = plt.plot(data['x'], data['y'], color=color, linewidth=linewidth)
                        plt.xlabel('$\\theta_%d$' % ix)
                        plt.ylabel('$\\theta_%d$' % iy)
                        if pct>0.7:

                            legend_handles.append(hh)
                            legend_labels.append(legend_dict[suffix][pct])

                    

    plt.legend(legend_handles, legend_labels, fontsize=10,
               bbox_to_anchor=(1.05, 1), loc=2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "truth_contours.png"))
    plt.close()


    plt.figsize = (30, 30)
    ix = 0
    iy = 1
    legend_handles = []
    legend_labels = []

    plt.subplot(2,2,1)
    for suffix, color in zip(('bayesianFullD', 'frequentistFullDrelative'),
                             ('b', 'm')):
                
        for pct in (0.68, 0.95):
                    
            file_name = '%s_%.2f_%d_%d_%s.txt' % \
            (control_root, pct, ix, iy, suffix)
                        
            data = np.genfromtxt(os.path.join(control_dir, file_name),
                                 dtype=dtype)

            ordered_x, ordered_y = order_x_y(data['x'], data['y'])
                
            label = legend_dict[suffix][pct]
            if pct>0.7:
                linewidth=2
            else:
                linewidth=0.5
                        
            hh, = plt.plot(ordered_x, ordered_y, color=color, linewidth=linewidth)
            plt.xlabel('$\\theta_%d$' % ix)
            plt.ylabel('$\\theta_%d$' % iy)
            if pct>0.7:
                legend_handles.append(hh)
                legend_labels.append(legend_dict[suffix][pct])

    plt.legend(legend_handles, legend_labels, fontsize=10,
               bbox_to_anchor=(1.05, 1), loc=2)


    legend_handles = []
    legend_labels = []

    plt.subplot(2,2,3)
    for suffix, color in zip(('bayesian2D', 'frequentist2D'),
                             ('r', 'g')):
                
        for pct in (0.68, 0.95):
                    
            file_name = '%s_%.2f_%d_%d_%s.txt' % \
            (control_root, pct, ix, iy, suffix)
                        
            data = np.genfromtxt(os.path.join(control_dir, file_name),
                                 dtype=dtype)

            ordered_x, ordered_y = order_x_y(data['x'], data['y'])
                
            label = legend_dict[suffix][pct]
            if pct>0.7:
                linewidth=2
            else:
                linewidth=0.5
                        
            hh, = plt.plot(ordered_x, ordered_y, color=color, linewidth=linewidth)
            plt.xlabel('$\\theta_%d$' % ix)
            plt.ylabel('$\\theta_%d$' % iy)
            if pct>0.7:
                legend_handles.append(hh)
                legend_labels.append(legend_dict[suffix][pct])



    plt.legend(legend_handles, legend_labels, fontsize=10,
               bbox_to_anchor=(1.05, 1), loc=2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "truth_contours_0_1.png"))
    plt.close()


    #plt.legend(legend_handles, legend_labels, fontsize=7, loc='lower right')
    #plt.tight_layout()
    #plt.savefig(os.path.join(output_dir, "truth_contours_%d.png" % iplot))
