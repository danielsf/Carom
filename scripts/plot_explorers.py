import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from analyzeCarom import scatter_from_carom
from analyzeCarom import scatter_from_multinest_projection

if __name__ == "__main__":

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    dalex_dir = os.path.join(physics_dir, "Carom", "output", "workspace")

    delta_chi = 21.03
    nlive = 50000

    multinest_file = os.path.join(multinest_dir,
                                  "gaussianJellyBean_d12_s99_n%d_t1.00e-03.txt" % nlive)

    dim_list = [(0,1), (6,9)]

    _data = None
    m_x_dict = {}
    m_y_dict = {}
    for dim in dim_list:
        m_x, m_y, _data = scatter_from_multinest_projection(multinest_file, 12, dim[0], dim[1],
                                                            data=_data)

        m_x_dict[dim] = m_x
        m_y_dict[dim] = m_y

    seed_list = [66, 694, 762, 1068, 6475, 626]

    data_dict = {}
    for seed in seed_list:
        data_dict[seed] = None

    dtype_list = []
    for ii in range(12):
        dtype_list.append(('p%d' % ii, float))
    dtype_list.append(('chisq', float))
    dtype_list.append(('j1', int))
    dtype_list.append(('j2', int))
    dtype_list.append(('log', int))

    dtype = np.dtype(dtype_list)

    for dim in dim_list:
        plt.figsize = (30, 30)
        for i_fig, seed in enumerate(seed_list):

            dalex_name = os.path.join(dalex_dir, "jellyBean_d12_s%d_output.sav" % seed)

            data = np.genfromtxt(dalex_name, dtype=dtype)
            ex_dex = np.where(data['log']==2)
            data = data[ex_dex]

            data_dict[seed] = data

            plt.subplot(3,2,i_fig+1)
            plt.scatter(m_x_dict[dim], m_y_dict[dim], color='k')
            plt.scatter(data['p%d' % dim[0]], data['p%d' % dim[1]], marker='x', color='r')
            plt.title('seed = %d' % (seed),
                      fontsize=10)

        plt.tight_layout()
        plt.savefig(os.path.join(physics_dir, 'Carom', 'figures',
                                 'dalex_explorers_%d_%d.png' % (dim[0],dim[1])))
        plt.close()
