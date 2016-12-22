import numpy as np

__all__ = ["load_multinest_data", "load_dalex_data",
           "scatter_from_multinest_projection",
           "scatter_from_multinest_marginalized",
           "scatter_from_carom", "make_histogram"]

def load_multinest_data(file_name, dim):
    dt_list = [('degen', np.float), ('chisq', np.float)]
    for ii in range(dim):
        dt_list.append(('x%d' % ii, np.float))
    dtype = dt_list
    return np.genfromtxt(file_name, dtype=dtype)


def load_dalex_data(data_name, dim):
    dt_list = []
    for ii in range(dim):
        dt_list.append(('x%d' % ii, np.float))

    dt_list.append(('chisq', np.float))
    dt_list.append(('junk1', int))
    dt_list.append(('junk2', int))
    dt_list.append(('junk3', int))
    dtype = np.dtype(dt_list)

    return np.genfromtxt(data_name, dtype=dtype)

def make_histogram(xx_in, dmag, cut_off, cumulative=True):
    xx = xx_in[np.where(xx_in<=cut_off+dmag)]
    #print xx.min(),xx.max()
    min_val=xx.min()-dmag
    i_xx = np.round((xx-min_val)/dmag).astype(int)
    unique_ixx, ct = np.unique(i_xx, return_counts=True)

    if cumulative:
        return unique_ixx*dmag+min_val, ct.astype(float)/float(len(xx_in))
    else:
        return unique_ixx*dmag+min_val, ct.astype(int)


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


def raw_bayes(data_x, data_y, density_in, prob=0.95):

    sorted_density = np.sort(density_in)
    total_sum = density_in.sum()
    sum_cutoff = 0.0
    for ii in range(len(sorted_density)-1, -1, -1):
        sum_cutoff += sorted_density[ii]
        if sum_cutoff >= prob*total_sum:
            cutoff = sorted_density[ii]
            break

    dexes = np.where(density_in>=cutoff)
    return data_x[dexes], data_y[dexes]


def get_scatter(data_x, data_y, x_norm, y_norm):

    x_out = np.ones(len(data_x))*data_x[-1]
    y_out = np.ones(len(data_y))*data_y[-1]

    tol = 0.0025

    i_x_min = np.argmin(data_x)
    i_x_max = np.argmax(data_x)
    i_y_min = np.argmin(data_y)
    i_y_max = np.argmax(data_y)

    actual_len = 1
    chosen_dexes = [len(data_x)-1]

    for dex in (i_x_min, i_x_max, i_y_min, i_y_max):
        if dex not in chosen_dexes:
            x_out[actual_len] = data_x[dex]
            y_out[actual_len] = data_y[dex]
            actual_len += 1
            chosen_dexes.append(dex)

    for ii in range(len(data_x)-1, -1, -1):
        dd_min = (np.power((data_x[ii]-x_out[:actual_len])/x_norm,2) +
                  np.power((data_y[ii]-y_out[:actual_len])/y_norm,2)).min()

        if dd_min > tol:
            x_out[actual_len] = data_x[ii]
            y_out[actual_len] = data_y[ii]
            actual_len += 1

    return x_out[:actual_len], y_out[:actual_len]


def scatter_from_multinest_projection(file_name, dim, ix, iy, data=None):

    if data is None:
        ref_data = load_multinest_data(file_name, dim)

    else:
        ref_data = data

    ref_x, ref_y = raw_bayes(ref_data['x%d' % ix],
                               ref_data['x%d' % iy],
                               ref_data['degen'])

    return ref_x, ref_y, ref_data

def scatter_from_multinest_marginalized(file_name, dim, ix, iy, data=None):

    if data is None:
        ref_data = load_multinest_data(file_name, dim)

    else:
        ref_data = data

    ref_x, ref_y = marginalize(ref_data['x%d' % ix],
                               ref_data['x%d' % iy],
                               ref_data['degen'])

    return ref_x, ref_y, ref_data


def scatter_from_carom(data_name, dim, ix, iy, delta_chi=None, target=None, data=None, limit=None):
    if data is None:
        data = load_dalex_data(data_name, dim)

    if limit is not None:
        data_cut = data[:limit]
    else:
        data_cut = data

    mindex = np.argmin(data_cut['chisq'])
    chisq_min = data_cut['chisq'][mindex]

    if delta_chi is not None:
        if target is not None:
            raise RuntimeWarning("You have specified both target and delta_chi; don't do that")
        target = chisq_min + delta_chi

    good_dexes = np.where(data_cut['chisq'] <= target)
    good_x = data_cut['x%d' % ix][good_dexes]
    good_y = data_cut['x%d' % iy][good_dexes]

    x_max = good_x.max()
    x_min = good_x.min()
    y_max = good_y.max()
    y_min = good_y.min()

    x_norm = x_max-x_min
    y_norm = y_max-y_min


    dd_arr = np.power((good_x-data_cut['x%d' % ix][mindex])/x_norm, 2) + \
             np.power((good_y-data_cut['x%d' % iy][mindex])/y_norm, 2)

    dd_sorted_dexes = np.argsort(dd_arr)

    x_grid, y_grid = get_scatter(good_x[dd_sorted_dexes],
                                 good_y[dd_sorted_dexes],
                                 x_norm, y_norm)

    return x_grid, y_grid, chisq_min, target, data
