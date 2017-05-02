import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import argparse

from analyzeDalex import scatter_from_dalex
from analyzeDalex import scatter_from_multinest_projection

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, default=None)
    parser.add_argument('--seed', type=int, nargs='+', default=None)
    parser.add_argument('--limit', type=int, default=600000)
    parser.add_argument('--input', type=str,
                        default='output/workspace')
    parser.add_argument('--x', type=int, nargs='+', default=None)
    args = parser.parse_args()

    if args.prefix is None:
        raise RuntimeError("must specify prefix")

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    dalex_dir = args.input

    delta_chi = 21.03
    nlive = 50000

    multinest_file = os.path.join(multinest_dir,
                                  "gaussianJellyBean_d12_s99_n%d_t1.00e-03.txt" % nlive)

    if args.x is None:
        dim_list = [(0,1), (6,9)]
    else:
        dim_list = []
        for ix in range(0,len(args.x),2):
            dim_list.append((args.x[ix], args.x[ix+1]))

    _data = None
    m_x_dict = {}
    m_y_dict = {}
    for dim in dim_list:
        m_x, m_y, _data = scatter_from_multinest_projection(multinest_file, 12, dim[0], dim[1],
                                                            data=_data)

        m_x_dict[dim] = m_x
        m_y_dict[dim] = m_y

    seed_list = args.seed
    
    data_dict = {}
    for seed in seed_list:
        data_dict[seed] = None
    
    for dim in dim_list:
        plt.figsize = (30, 30)
        for i_fig, seed in enumerate(seed_list):
                
            dalex_name = os.path.join(dalex_dir, "jellyBean_d12_s%d_output.sav" % seed)
            
            (d_x, d_y, d_min, d_target,
             data) = scatter_from_dalex(dalex_name, 12, dim[0], dim[1],
                                        delta_chi=delta_chi, data=data_dict[seed],
                                        limit=args.limit)
                
            data_dict[seed] = data
                
            plt.subplot(3,2,i_fig+1)
            plt.scatter(m_x_dict[dim], m_y_dict[dim], color='k')
            plt.scatter(d_x, d_y, marker='x', color='r')
            plt.title('$\chi^2_{min} = %.2f$; n_calls = %.2e\nseed %d' % (d_min, args.limit, seed),
                      fontsize=15)
            
        plt.tight_layout()
        plt.savefig('%s_%d_%d_%d.png' % (args.prefix,args.limit,dim[0],dim[1]))
        plt.close()
