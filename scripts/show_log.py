import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from analyzeCarom import scatter_from_carom, make_histogram

if __name__ == "__main__":

    limit = 250000
    ix = 6
    iy = 9
    dim = 12
    delta =21.03

    type_dict = {}
    type_dict['init'] = 0
    type_dict['refine'] = 1
    type_dict['explore'] = 2
    type_dict['tendril'] = 3

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    carom_dir = os.path.join(physics_dir, "Carom")
    data_file = os.path.join(carom_dir, "output", "draft_161215",
                             "jellyBean_log_output.sav")

    (final_x, final_y, chisq_min,
     chisq_target, data) = scatter_from_carom(data_file, dim, ix, iy,
                                              delta_chi=delta, limit=limit)

    dtype_list = []
    for ii in range(dim):
        dtype_list.append(('p%d' % ii, float))
    dtype_list.append(('chisq', float))
    dtype_list.append(('mu', int))
    dtype_list.append(('sig', int))
    dtype_list.append(('log', int))

    dtype = np.dtype(dtype_list)
    raw_data = np.genfromtxt(data_file, dtype=dtype)
    raw_data = raw_data[:limit]
    assert len(raw_data) == limit

    for type_key in type_dict:
        scatter_name = os.path.join(fig_dir, "time_evolution_%s.png" % type_key)
        i_type =type_dict[type_key]
        
        valid_rows = np.where(raw_data['log'] == i_type)
        xx_master = raw_data['p%d' % ix][valid_rows]
        yy_master = raw_data['p%d' % iy][valid_rows]
        chisq_master = raw_data['chisq'][valid_rows]
        print type_key,' chimax ',chisq_master.max(),chisq_master.min()
        
        _xmax = max(xx_master.min(), xx_master.max(), final_x.max(), final_x.min())
        _xmin = min(xx_master.min(), xx_master.max(), final_x.max(), final_x.min())

        _ymax = max(yy_master.min(), yy_master.max(), final_y.max(), final_y.min())
        _ymin = min(yy_master.min(), yy_master.max(), final_y.max(), final_y.min())
        
        _dx = 0.01*(_xmax-_xmin)
        _dy = 0.01*(_ymax-_ymin)
        
        plt.figsize = (30, 30)
        for i_fig in range(1,5):
            plt.subplot(2, 2, i_fig)
            cut_off = i_fig*len(xx_master)/4-1
            xx = xx_master[:cut_off]
            yy = yy_master[:cut_off]
            print cut_off, len(valid_rows[0])
            tt = valid_rows[0][cut_off]
            
            plt.scatter(final_x, final_y, color='y', s=10)
            plt.scatter(xx, yy, color='r', marker='x', s=5)
            plt.xlim((_xmin, _xmax))
            plt.ylim((_ymin, _ymax))
            plt.text(_xmin+_dx, _ymin+_dy,
                     '%d total $\chi^2$ calls' % (tt),
                     fontsize=10)
        
        
        plt.savefig(scatter_name)
        plt.close()
