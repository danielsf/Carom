import os
import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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
    ix = 6
    iy = 9
    seeds = [234, 786, 932] #, 99, 66, 125, 6475]

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

    total_post = ref_data['degen'].sum()
    sum_post = 0.0
    for ii in range(len(ref_data)-1, -1, -1):
        sum_post += ref_data['degen'][ii]
        cutoff = ref_data['degen'][ii]
        if sum_post >= 0.95*total_post:
            break

    good_dexes = np.where(ref_data['degen']>cutoff)
    raw_ref_x = ref_data['x%d' % ix][good_dexes]
    raw_ref_y = ref_data['x%d' % iy][good_dexes]

    x_max = raw_ref_x.max()
    x_min = raw_ref_x.min()
    y_max = raw_ref_y.max()
    y_min = raw_ref_y.min()

    x_norm = x_max-x_min
    y_norm = y_max-y_min

    dd_arr = np.power((raw_ref_x - 0.5*(x_max+x_min))/x_norm,2) + \
             np.power((raw_ref_y - 0.5*(y_max+y_min))/y_norm,2)

    sorted_dexes = np.argsort(dd_arr)

    ref_x, ref_y = get_scatter(raw_ref_x[sorted_dexes],
                               raw_ref_y[sorted_dexes],
                               x_norm, y_norm)

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
            out_name = os.path.join(data_dir, "jb_matrix_%d.eps" % ifig)


        n_row = ct/col_max
        n_col = ct%row_max
        ax = ax_arr[n_row][n_col]
        data_name = os.path.join(data_dir,'jellyBean_d%d_s%d_output.sav'
                                 % (dim, ss))

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
        ax.text(x_max-20, y_max-5,
               '$\chi^2_{min}$=%.2f\npts %d' % (chisq_min, len(data['x0'])),
               fontdict={'fontsize':10})


        ct += 1
        if ct == row_max*col_max or ss==seeds[-1]:
            fig.savefig(out_name)
            fig = None


