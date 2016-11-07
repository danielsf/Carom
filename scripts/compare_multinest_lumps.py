import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os

from analyzeCarom import scatter_from_multinest_projection

import sys

if __name__ == "__main__":

    dim = (int(sys.argv[1]), int(sys.argv[2]))
    n_list = []
    for ii in range(3, len(sys.argv)):
        n_list.append(int(sys.argv[ii]))

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")

    multinest_carom_name = [os.path.join(multinest_dir,
                                        "nonGaussianLump_d12_s99_n%d_t1.00e-03.txt" % ii)
                            for ii in n_list]

    plt.figsize = (30,30)
    xmin = 1.0e30
    xmax = -1.0e30
    ymin = 1.0e30
    ymax = -1.0e30
    full_dim = 12
    color_list = ['k', 'r', 'b', 'c', 'g']
    hh_list = []
    tag_list = []
    for i_fig in range(len(multinest_carom_name)):
        name = multinest_carom_name[i_fig]
        xx, yy, data = scatter_from_multinest_projection(name, full_dim, dim[0], dim[1])
        hh = plt.scatter(xx, yy, color=color_list[i_fig])
        hh_list.append(hh)
        tag_list.append('%d live points' % n_list[i_fig])
        print xx.min(),xx.max(),len(xx)
        print yy.min(),yy.max(),len(yy)
        for vv in (xx.min(), xx.max()):
            if vv<xmin:
                xmin=vv
            if vv>xmax:
                xmax=vv
        for vv in (yy.min(), yy.max()):
            if vv<ymin:
                ymin=vv
            if vv>ymax:
                ymax=vv
    
    dx=xmax-xmin
    dy=ymax-ymin
    print xmin,xmax
    print ymin,ymax
    plt.xlim((xmin-0.05*dx, xmax+0.05*dx))
    plt.ylim((ymin-0.05*dy, ymax+0.05*dy))
    plt.legend(hh_list, tag_list, fontsize=10, loc=0)
    plt.xlabel('$\\theta_%d$' % dim[0], fontsize=20)
    plt.ylabel('$\\theta_%d$' % dim[1], fontsize=20)
    plt.savefig(os.path.join(fig_dir,'lump_mult_compare_%d_%d.png' % (dim[0], dim[1])))
    plt.close()
    
    
