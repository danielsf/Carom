from __future__ import with_statement
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import os

if __name__ == "__main__":

    scratch_dir = os.path.join('/Users', 'danielsf', 'physics', 'Carom')
    scratch_dir = os.path.join(scratch_dir, 'controls', 'scratch')

    ix_list = (0, 0, 0, 1, 1, 2)
    iy_list = (1, 2, 3, 2, 3, 3)

    data = {}

    for ix, iy in zip(ix_list, iy_list):
        name = 'integrable_extended_0.95_%d_%d_heatmap.txt' % (ix, iy)
        input_name = os.path.join(scratch_dir, name)
        data['%d_%d' % (ix, iy)] = np.genfromtxt(input_name).transpose()


    plt.figure(figsize=(30, 30))
    for ifig, (ix, iy) in enumerate(zip(ix_list, iy_list)):

        tag = '%d_%d' % (ix, iy)
        plt.subplot(3, 2, ifig+1)
        indices = np.where(data[tag][2]<50.0)

        xx = data[tag][0][indices]
        yy = data[tag][1][indices]
        cc = data[tag][2][indices]

        x_min = xx.min()
        x_max = xx.max()
        y_min = yy.min()
        y_max = yy.max()

        xx_sorted = np.sort(xx)
        dx_raw = np.unique(np.diff(xx_sorted))
        dx = dx_raw[np.where(dx_raw>1.0e-10)].min()

        yy_sorted = np.sort(yy)
        dy_raw = np.unique(np.diff(yy_sorted))
        dy = dy_raw[np.where(dy_raw>1.0e-10)].min()

        xx_grid, yy_grid = np.mgrid[slice(x_min, x_max+dx, dx),
                                    slice(y_min, y_max+dy, dy)]

        cc_grid = np.ones(xx_grid.shape)*(-99.0)

        for ii in range(len(xx)):
            ix = int((xx[ii]-x_min)/dx)
            iy = int((yy[ii]-y_min)/dy)
            cc_grid[ix][iy] = cc[ii]

        cc_masked = np.ma.masked_values(cc_grid, -99.0)


        plt.pcolormesh(xx_grid, yy_grid, cc_masked,
                       cmap=plt.cm.gist_ncar,
                       edgecolor='',
                       norm=LogNorm(vmin=cc.min(), vmax=cc.max()))

        ticks=np.linspace(0.0,30.0,5.0,endpoint=True)
        #plt.contour(data[tag][0], data[tag][1], data[tag][2],ticks)
        xmax=data[tag][0][indices].max()
        xmin=data[tag][0][indices].min()
        ymax=data[tag][1][indices].max()
        ymin=data[tag][1][indices].min()
        dx=(xmax-xmin)*0.1
        dy=(ymax-ymin)*0.1
        plt.text(xmax-2*dx, ymax-1*dy, tag, fontsize=30)
        plt.colorbar()

    plt.savefig(os.path.join(scratch_dir, 'heat_map.eps'))
