from __future__ import with_statement
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

if __name__ == "__main__":

    scratch_dir = os.path.join('/Users', 'danielsf', 'physics', 'Carom')
    scratch_dir = os.path.join(scratch_dir, 'controls', 'scratch')

    ix_list = (0, 0, 0, 1, 1, 2)
    iy_list = (1, 2, 3, 2, 3, 3)

    data = {}
    
    for ix, iy in zip(ix_list, iy_list):
        name = 'integrable_gross_0.95_%d_%d_heatmap.txt' % (ix, iy)
        input_name = os.path.join(scratch_dir, name)
        data['%d_%d' % (ix, iy)] = np.genfromtxt(input_name).transpose()
        

    plt.figure(figsize=(30, 30))
    for ifig, (ix, iy) in enumerate(zip(ix_list, iy_list)):
        tag = '%d_%d' % (ix, iy)
        plt.subplot(3, 2, ifig+1)
        plt.scatter(data[tag][0], data[tag][1],
                    c=data[tag][2],
                    #c=np.exp(-0.5*data[tag][2]/data[tag][2].max()),
                    cmap=plt.cm.gist_ncar,
                    edgecolor='')
        xmax=data[tag][0].max()
        xmin=data[tag][0].min()
        ymax=data[tag][1].max()
        ymin=data[tag][1].min()
        dx=(xmax-xmin)*0.1
        dy=(ymax-ymin)*0.1
        plt.text(xmax-2*dx, ymax-1*dy, tag, fontsize=30)
        plt.colorbar()

    plt.savefig(os.path.join(scratch_dir, 'heat_map.eps'))
