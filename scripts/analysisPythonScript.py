from __future__ import with_statement
import numpy as np

def get_histogram(chisq, dchi, nsteps):
    chisq_min = chisq.min()
    x_arr = np.arange(chisq_min+0.5*dchi, chisq_min+(nsteps+1)*dchi, dchi)
    cts = np.zeros(len(x_arr))
    for ix, xx in enum(x_arr):
        cts[ix] = len(np.where(chisq<xx+0.5*dchi)[0])
        if ix>0:
            cts[ix]-=cts[ix-1]

    return x_arr, cts


def get_good_pts(data, target=None, delta_chisq=None):
    if target is None and delta_chisq is None:
        raise RuntimeError("must specify either target or delta_chisq in get_scatter")


    i_chi = data.shape[1]-1

    if target is not None:
        _target = target
    else:
        _target = data.transpose()[i_chi].min()+delta_chisq


    good_pts = data[np.where(data.transpose()[i_chi]<=_target)[0]]
    return good_pts


def get_scatter(data, ix, iy, target=None, delta_chisq=None, ddsq_threshold=0.001):

    good_pts =get_good_pts(data, target=target, delta_chisq=delta_chisq)

    i_chi = data.shape[1]-1

    temp = good_pts.transpose()
    x_norm = temp[ix].max()-temp[ix].min()
    y_norm = temp[iy].max()-temp[iy].min()

    min_dex = np.argmin(temp[i_chi])
    print temp[i_chi][min_dex]
    print len(temp[i_chi])
    ct_kept = 1
    pts_kept = np.ones((2,len(data)))*1000000.0
    pts_kept[0][0] = good_pts[min_dex][ix]
    pts_kept[1][0] = good_pts[min_dex][iy]

    chisq_kept = [good_pts[min_dex][i_chi]]
    print 'data shape ',data.shape,i_chi
    print good_pts[min_dex]
    print 'min_dex ',min_dex

    print 'looping ',x_norm, y_norm
    for pt in good_pts:
        print_it = False

        if not print_it:
            dd = np.power((pt[ix]-pts_kept[0][:ct_kept])/x_norm,2)+np.power((pt[iy]-pts_kept[1][:ct_kept])/y_norm,2)
            if dd.min()>ddsq_threshold:
                print_it = True

        if print_it:
            pts_kept[0][ct_kept] = pt[ix]
            pts_kept[1][ct_kept] = pt[iy]
            chisq_kept.append(pt[i_chi])
            ct_kept += 1

    if len(chisq_kept)!=ct_kept:
        raise RuntimeError("ct_kept %d chisq %d\n" %(ct_kept, len(chisq_kept)))
    return pts_kept.transpose()[:ct_kept], chisq_kept


import os

if __name__ == "__main__":

    dim = 10

    input_dir = os.path.join('/Users','danielsf','physics','Carom','output','scratch')

    input_file = os.path.join(input_dir,'test151023','jellyBean_d10_s99_output.sav')

    full_data = np.genfromtxt(input_file)
    print full_data[dim].min()
    useful_data = full_data.transpose()[:dim+1].transpose()
    print useful_data[0]
    print full_data.shape
    print useful_data.shape
    print full_data.view(np.ndarray)

    good_pts, chisq = get_scatter(useful_data, 0, 1, delta_chisq=18.3)
    print 'final ',good_pts.shape
    with open(os.path.join(input_dir,'trial_plot.sav'), 'w') as output:
        for pt, cc in zip(good_pts, chisq):
            output.write('%e %e %e\n' % (pt[0], pt[1], cc))
