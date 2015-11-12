from __future__ import with_statement
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

__all__ = ["doAnalysis"]

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
        raise RuntimeError("must specify either target or delta_chisq in get_good_pts")

    i_chi = data.shape[1]-1

    if target is not None:
        _target = target
    else:
        _target = data.transpose()[i_chi].min()+delta_chisq


    good_pts = data[np.where(data.transpose()[i_chi]<=_target)[0]]

    return good_pts, _target


def _get_scatter(good_pts, ix, iy, ddsq_threshold=0.001):

    i_chi = good_pts.shape[1]-1

    temp = good_pts.transpose()
    x_norm = temp[ix].max()-temp[ix].min()
    y_norm = temp[iy].max()-temp[iy].min()

    min_dex = np.argmin(temp[i_chi])
    ct_kept = 1
    pts_kept = np.ones((2,len(good_pts)))*1000000.0
    pts_kept[0][0] = good_pts[min_dex][ix]
    pts_kept[1][0] = good_pts[min_dex][iy]

    chisq_kept = [good_pts[min_dex][i_chi]]

    center_dd = np.power((temp[ix]-pts_kept[0][0])/x_norm,2)+np.power((temp[iy]-pts_kept[1][0])/y_norm,2)

    for ipt in range(len(good_pts)):
        print_it = False

        i_chosen = np.argmax(center_dd)

        pt = good_pts[i_chosen]

        if not print_it:
            dd = np.power((pt[ix]-pts_kept[0][:ct_kept])/x_norm,2)+np.power((pt[iy]-pts_kept[1][:ct_kept])/y_norm,2)
            if dd.min()>ddsq_threshold:
                print_it = True

        if print_it:
            pts_kept[0][ct_kept] = pt[ix]
            pts_kept[1][ct_kept] = pt[iy]
            chisq_kept.append(pt[i_chi])
            ct_kept += 1

        center_dd[i_chosen] = -1000.0

    if len(chisq_kept)!=ct_kept:
        raise RuntimeError("ct_kept %d chisq %d\n" %(ct_kept, len(chisq_kept)))
    return pts_kept.transpose()[:ct_kept], chisq_kept


def get_scatter(data, ix, iy, target=None, delta_chisq=None, ddsq_threshold=0.001):

    good_pts, _target =get_good_pts(data, target=target, delta_chisq=delta_chisq)
    out_tuple = _get_scatter(good_pts, ix, iy, ddsq_threshold=ddsq_threshold)
    return out_tuple[0], out_tuple[1], _target

import os


def doAnalysis(dim, delta_chisq, ix_list, iy_list, ct_list, input_file, control_names, output_dir, scatter_control=False):

    plot_rows = len(ct_list)/3
    if 3*plot_rows<len(ct_list):
        plot_rows += 1


    full_data = np.genfromtxt(input_file)
    useful_data = full_data.transpose()[:dim+1].transpose()

    log_suffixes = {'ricochet':'_ricochet_log.txt',
                    'simplex':'_simplex_log.txt',
                    'dchi_simplex':'_dchi_simplex_log.txt',
                    'mcmc':'_mcmc_log.txt',
                    'swarm':'_swarm_log.txt',
                    'compass':'_compass_log.txt'}


    log_data = {}
    for log_name in log_suffixes:
        log_file = input_file + log_suffixes[log_name]
        _log_data = np.genfromtxt(log_file)
        log_data[log_name] = _log_data


    for ix, iy in zip(ix_list, iy_list):

        control_data = np.genfromtxt(control_names['%d_%d' % (ix, iy)]).transpose()

        xmax = control_data[0].max()
        xmin = control_data[0].min()
        ymax = control_data[1].max()
        ymin = control_data[1].min()

        dx = (xmax-xmin)/9.0
        xmin -= dx
        xmax += dx

        dy = (ymax-ymin)/9.0
        ymin -= dy
        ymax += dy

        xticks = np.arange(xmin,xmax+dx,dx)
        xformat = ['%.2e' % xticks[ii] if ii%3==0 else '' for ii in range(len(xticks))]

        yticks = np.arange(ymin,ymax+dx,dy)
        yformat = ['%.2e' % yticks[ii] if ii%3==0 else '' for ii in range(len(yticks))]

        plt.figure(figsize=(30,30))
        for ict, ct in enumerate(ct_list):
            all_pts = useful_data[:ct]
            good_pts, chisq, target = get_scatter(useful_data[:ct], ix, iy, delta_chisq=delta_chisq)

            good_pts = good_pts.transpose()

            plt.subplot(plot_rows, 3, ict+1)
            if not scatter_control:
                plt.plot(control_data[0], control_data[1], linewidth=2)
            else:
                plt.scatter(control_data[0], control_data[1], color='b', marker='x', s=40)
            plt.scatter(good_pts[0], good_pts[1], color='r', s=40)

            plt.xlabel('$\\theta_%d$' % ix, fontsize=20)
            plt.ylabel('$\\theta_%d$' % iy, fontsize=20)
            plt.xticks(xticks, xformat, fontsize=20)
            plt.yticks(yticks, yformat, fontsize=20)

            title = 'target $\chi^2 =$ %.2f\npoints %d' % (target, ct)

            plt.text(xmax-0.6*(xmax-xmin), ymax-0.2*(ymax-ymin), title, fontsize=30)

        file_name = os.path.join(output_dir, 'full_%d_%d.eps' % (ix, iy))

        plt.savefig(file_name)
        plt.close()

        for log_name in ('ricochet', 'swarm', 'dchi_simplex', 'mcmc', 'compass'):
            temp = log_data[log_name].transpose()
            data_dexes = temp[dim+1]
            plt.figure(figsize=(30,30))
            for ict, ct in enumerate(ct_list):
                chosen_dexes = np.where(data_dexes<=ct)[0]
                data = log_data[log_name][chosen_dexes].transpose()
                plt.subplot(plot_rows, 3, ict+1)
                if not scatter_control:
                    plt.plot(control_data[0], control_data[1], linewidth=2)
                else:
                    plt.scatter(control_data[0], control_data[1], color='b', marker='x', s=40)
                plt.scatter(data[ix], data[iy], color='r', s=40)
                plt.xlabel('$\\theta_%d$' % ix, fontsize=20)
                plt.ylabel('$\\theta_%d$' % iy, fontsize=20)
                plt.xticks(xticks, xformat, fontsize=20)
                plt.yticks(yticks, yformat, fontsize=20)

                title = 'points %d' % (ct)
                plt.text(xmax-0.6*(xmax-xmin), ymax-0.2*(ymax-ymin), title, fontsize=30)

            file_name = os.path.join(output_dir, '%s_%d_%d.eps' % (log_name, ix, iy))
            plt.savefig(file_name)
            plt.close()