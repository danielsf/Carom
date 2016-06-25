import os
import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def marginalize(data_x, data_y, density_in, prob=0.95):

    i_dx = 100

    xmax = data_x.max()
    xmin = data_x.min()
    dx = (xmax-xmin)/float(i_dx)

    ymax = data_y.max()
    ymin = data_y.min()
    dy = (ymax-ymin)/float(i_dx)

    i_x = np.array(range((i_dx+1)*(i_dx+1)))

    x_out = np.array([xmin + (ii%i_dx)*dx for ii in i_x])
    y_out = np.array([ymin + (ii/i_dx)*dy for ii in i_x])

    density_out = np.zeros(len(x_out))

    for xx, yy, dd in zip(data_x, data_y, density_in):
        ix = int(round((xx-xmin)/dx))
        iy = int(round((yy-ymin)/dy))
        dex = ix + iy*i_dx
        density_out[dex] += dd

    sorted_density = np.sort(density_out)
    total_sum = density_out.sum()
    sum_cutoff = 0.0
    for ii in range(len(sorted_density)-1, -1, -1):
        sum_cutoff += sorted_density[ii]
        if sum_cutoff >= prob*total_sum:
            cutoff = sorted_density[ii]
            break

    dexes = np.where(density_out>=cutoff)
    return x_out[dexes], y_out[dexes]


def get_scatter(data_x, data_y, x_norm, y_norm):

    x_out = np.ones(len(data_x))*data_x[-1]
    y_out = np.ones(len(data_y))*data_y[-1]

    tol = 0.0025

    actual_len = 1
    for ii in range(len(data_x)-1, -1, -1):
        dd_min = (np.power((data_x[ii]-x_out[:actual_len])/x_norm,2) +
                  np.power((data_y[ii]-y_out[:actual_len])/y_norm,2)).min()

        if dd_min > tol:
            x_out[actual_len] = data_x[ii]
            y_out[actual_len] = data_y[ii]
            actual_len += 1

    return x_out[:actual_len], y_out[:actual_len]



if __name__ == "__main__":

    dim = 12
    delta_chi = 21.0
    ix = 0
    iy = 1
    seeds = [234, 786, 932, 99, 66, 125, 6475]
    fig_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(fig_dir, "Carom", "figures")
    fig_dir = os.path.join(fig_dir, "dd")

    row_max = 2
    col_max = 2

    data_dir = os.path.join("/Users", "danielsf", "physics")
    data_dir = os.path.join(data_dir, "Carom", "output", "scratch")

    ref_dir = os.path.join("/Users", "danielsf", "physics")
    ref_dir = os.path.join(ref_dir,"Multinest_v3.9", "chains")

    ref_file = os.path.join(ref_dir, "gaussianJellyBean_d12_s99_n300.txt")

    dt_list = [('degen', np.float), ('chisq', np.float)]
    for ii in range(dim):
        dt_list.append(('x%d' % ii, np.float))
    dtype = dt_list
    ref_data = np.genfromtxt(ref_file, dtype=dtype)

    ref_x, ref_y = marginalize(ref_data['x%d' % ix],
                               ref_data['x%d' % iy],
                               ref_data['degen'])

    dt_list = []
    for ii in range(dim):
        dt_list.append(('x%d' % ii, np.float))

    dt_list.append(('chisq', np.float))
    dt_list.append(('junk1', int))
    dt_list.append(('junk2', int))
    dt_list.append(('junk3', int))
    dtype = np.dtype(dt_list)
    fig = None
    ifig = -1
    t_start = time.clock()
    for ss in seeds:
        print "seed ",ss,time.clock()-t_start
        if fig is None:
            ct = 0
            ifig += 1
            fig, ax_arr = plt.subplots(nrows=row_max, ncols=col_max)
            fig.figsize=(30,30)
            out_name = os.path.join(fig_dir, "jb_dd_matrix_%d_x%d_%d.eps" %
                                    (ifig, ix, iy))


        n_row = ct/col_max
        n_col = ct%row_max
        ax = ax_arr[n_row][n_col]
        data_name = os.path.join(data_dir,'jellyBean_d%d_s%d_dd_output.sav'
                                 % (dim, ss))

        if os.path.exists(data_name):

            data = np.genfromtxt(data_name, dtype=dtype)
            mindex = np.argmin(data['chisq'])
            chisq_min = data['chisq'][mindex]
            target = chisq_min + delta_chi

            good_dexes = np.where(data['chisq'] <= target)
            good_x = data['x%d' % ix][good_dexes]
            good_y = data['x%d' % iy][good_dexes]


            x_max = good_x.max()
            x_min = good_x.min()
            y_max = good_y.max()
            y_min = good_y.min()

            x_norm = x_max-x_min
            y_norm = y_max-y_min


            dd_arr = np.power((good_x-data['x%d' % ix][mindex])/x_norm, 2) + \
                     np.power((good_y-data['x%d' % iy][mindex])/y_norm, 2)

            dd_sorted_dexes = np.argsort(-1.0*dd_arr)

            x_grid, y_grid = get_scatter(good_x[dd_sorted_dexes],
                                         good_y[dd_sorted_dexes],
                                         x_norm, y_norm)

            ax.scatter(ref_x, ref_y, color = 'k', s=5)
            ax.scatter(x_grid, y_grid, color='r', s=5)
            ax.set_title('seed: %d' % ss, fontdict={'fontsize':10})
            ax.text(x_min+0.75*(x_max-x_min), y_min+0.75*(y_max-y_min),
                   '$\chi^2_{min}$=%.2f\npts %d' % (chisq_min, len(data['x0'])),
                   fontdict={'fontsize':10})


        ct += 1
        if ct == row_max*col_max or ss==seeds[-1]:
            fig.savefig(out_name)
            fig = None


