from __future__ import with_statement
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

import scipy.spatial as spatial

from analyzeCarom import scatter_from_multinest_projection, scatter_from_carom
from analyzeCarom import scatter_from_multinest_marginalized

def sales_metric(pts):
    ans = np.sqrt(np.power(pts[:-1]-pts[1:],2).sum(axis=1)).sum()
    return ans

def replace(i1, i2, pts, trial_pts):
    if i1<i2:
        imin = i1
        imax = i2
    else:
        imin = i2
        imax = i1

    if imax-imin == 1:
        trial_pts[imin][0]=pts[imax][0]
        trial_pts[imin][1]=pts[imax][1]
        trial_pts[imax][0]=pts[imin][0]
        trial_pts[imax][1]=pts[imin][1]
    else:
        for ii in range(imin+1,imax+1):
            trial_pts[ii-1][0] = trial_pts[ii][0]
            trial_pts[ii-1][1] = trial_pts[ii][1]
        trial_pts[imax][0] = pts[imin][0]
        trial_pts[imax][1] = pts[imin][1]



def refine_salesman(rng, pts, iterations):

    orig_shape = pts.shape

    tree = spatial.KDTree(pts)
    n_nn = len(pts)/2

    print 'pts shape ',pts.shape

    d_current = sales_metric(pts)
    print 'start with ',d_current

    i_last = 0
    bad_pts = None
    for i_iteration in range(iterations):
        if bad_pts is None:
            dist, dex = tree.query(pts, 3)
            dex = dex.transpose()
            bad_pts = np.where(np.logical_and(np.abs(dex[1]-dex[0])>1, np.abs(dex[2]-dex[0])>1))[0]

        i1 = bad_pts[rng.random_integers(0,len(bad_pts)-1,1)[0]]

        dist, dex = tree.query(pts[i1], n_nn)
        has_chosen = []
        has_changed = False
        i_attempt = 0
        while not has_changed and i_attempt<10:
            i_attempt += 1
            i2 = i1
            while i2==i1 or i2 in has_chosen:
                i2 = dex[rng.random_integers(0,n_nn-1,1)[0]]
            has_chosen.append(i2)
            trial_pts = np.copy(pts)
            replace(i1, i2, pts, trial_pts)

            d_trial = sales_metric(trial_pts)
            if d_trial<d_current:
                d_current = d_trial
                pts = np.copy(trial_pts)
                tree = spatial.KDTree(pts)
                i_last = i_iteration
                bad_pts = None
                has_changed = True

        if i_iteration-i_last>1000:
            break


        print d_current
    print pts.shape,orig_shape
    return pts

def boundary_from_multinest_projection(file_name, dim, ix, iy, data=None):
    t_start = time.time()

    m_x, m_y, data = scatter_from_multinest_projection(file_name, dim, ix, iy, data=data)

    xmin = m_x.min()
    xmax = m_x.max()
    ymin = m_y.min()
    ymax = m_y.max()
    dx = 0.05*(xmax-xmin)
    dy = 0.05*(ymax-ymin)

    m_d_x = np.round(m_x/dx)
    m_d_y = np.round(m_y/dy)
    tree_data = np.array([m_d_x, m_d_y]).transpose()
    dex_raw = np.ascontiguousarray(tree_data).view(np.dtype((np.void, tree_data.dtype.itemsize*tree_data.shape[1])))
    _, unique_rows = np.unique(dex_raw, return_index=True)
    tree_data = tree_data[unique_rows]
    tree =spatial.KDTree(tree_data)
    boundary_pts = []
    for pt in tree_data:
        dist, dex = tree.query(pt, 5)
        if dist[4]>1.001:
            boundary_pts.append([pt[0], pt[1]])

    print "scatter %d boundary %d" % (len(tree_data), len(boundary_pts))

    boundary_pts_raw = np.array(boundary_pts).transpose()
    return boundary_pts_raw[0]*dx, boundary_pts_raw[1]*dy, data


    boundary_pts = []
    boundary_pts.append(np.copy(boundary_pts_raw[0]))
    boundary_pts_raw = np.delete(boundary_pts_raw, 0, 0)
    while len(boundary_pts_raw)>0:
        last_dex = len(boundary_pts)-1
        dex = None
        if len(boundary_pts_raw)==1:
            dex = 0
        if dex is None:
            possibilities = []
            dd_possibilities = []
            dd_last = np.sqrt(np.power(boundary_pts_raw-boundary_pts[last_dex],2).sum(axis=1))
            tree = spatial.KDTree(boundary_pts_raw)
            dist_arr, dex_arr = tree.query(boundary_pts_raw, 2)
            i_poss = np.where(dd_last<dist_arr.transpose()[1])
            dd_possibilities = dd_last[i_poss]
            possibilities = i_poss[0]

            if len(possibilities)>0:
                dex=possibilities[np.argsort(dd_possibilities)[0]]

        if dex is None:
            distance = np.sqrt(np.power(boundary_pts[last_dex]-boundary_pts_raw,2).sum(axis=1))
            dex = np.argmin(distance)

        pt = np.copy(boundary_pts_raw[dex])
        boundary_pts.append(pt)
        boundary_pts_raw = np.delete(boundary_pts_raw, dex, 0)

    un_refined = np.array(boundary_pts).transpose()
    print "boundary took ",time.time()-t_start
    return un_refined[0]*dx, un_refined[1]*dy, data

    """
    #boundary_pts = refine_salesman(rng, np.array(boundary_pts), 50000)

    #print boundary_pts[0]
    #print boundary_pts[1]

    for pt in un_refined.transpose():
        dd = np.power(pt-boundary_pts,2).sum(axis=1).min()
        try:
            assert dd<1.0e-30
        except:
            dex = np.argmin(np.power(pt-boundary_pts,2).sum(axis=1))
            print pt
            print boundary_pts[dex]
            print dd
            raise

    boundary_pts = boundary_pts.transpose()



    plt.figsize=(30,30)
    plt.scatter(m_x, m_y, color='k')
    #plt.plot(boundary_pts[0]*dx, boundary_pts[1]*dy, color='b', linestyle='-')
    plt.plot(un_refined[0]*dx, un_refined[1]*dy, color='r', linestyle='-')
    plt.savefig('junk.png')
    return un_refined[0], un_refined[1], data
    """


if __name__ == "__main__":

    seed_list = [66, 694, 762, 1068, 6475]
    limit_list = [400000, 300000, 400000, 400000, 400000]
    title_list = ['(a)', '(b)', '(c)', '(d)', '(e)']
    color_list = ['r', 'y', 'g', 'c']

    delta_chi = 21.03
    full_dim = 12

    physics_dir = os.path.join("/Users", "danielsf", "physics")
    fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")
    dalex_dir = os.path.join(physics_dir, "Carom")

    control_dir = os.path.join(dalex_dir, "controls", "draft_160907")
    multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
    data_dir = os.path.join(dalex_dir, "output", "draft_161117")

    dim_list = [(6,9), (0,1)]

    d_data_dict = {}
    m_data_dict = {}
    for seed in seed_list:
        d_data_dict[seed] = None
        m_data_dict[seed] = None


    for dim in dim_list:
        plt.figsize = (30,30)

        text_list = []
        xmax = None
        xmin = None
        ymax = None
        ymin = None

        for i_seed, (title, seed, limit) in enumerate(zip(title_list, seed_list, limit_list)):
            t_start = time.time()
            dalex_name = "jellyBean_d12_s%d_output.sav" % seed
            multinest_name = "gaussianJellyBean_d12_s%d_n300_t1.00e-03.txt" % seed

            scatter_name = "gaussianJellyBean_d12_s%d_n300_t1.00e-03_carom.sav" % seed
            with open(os.path.join(multinest_dir, scatter_name), 'r') as input_file:
                lines = input_file.readlines()
                n_mult = len(lines)

            plt.subplot(3,2,i_seed+1)
            plt.title(title, fontsize=7)

            d_h_list = []
            d_label_list = []

            (d_x, d_y, chisq_min, target,
             d_data) = scatter_from_carom(os.path.join(data_dir, dalex_name),
                                          full_dim, dim[0], dim[1], delta_chi=delta_chi,
                                          data=d_data_dict[seed],
                                          limit=limit)

            d_data_dict[seed] = d_data

            (d_x_forced, d_y_forced, chisq_min_forced, target_forced,
             d_data_forced) = scatter_from_carom(os.path.join(data_dir, dalex_name),
                                          full_dim, dim[0], dim[1], target=116.03,
                                          data=d_data_dict[seed],
                                          limit=limit)


            if xmax is None or d_x.max()>xmax:
                xmax=d_x.max()
            if xmin is None or d_x.min()<xmin:
                xmin=d_x.min()
            if ymax is None or d_y.max()>ymax:
                ymax=d_y.max()
            if ymin is None or d_y.min()<ymin:
                ymin=d_y.min()

            m_x, m_y, m_data = scatter_from_multinest_projection(
                                 os.path.join(multinest_dir, multinest_name),
                                 full_dim, dim[0], dim[1],
                                 data=m_data_dict[seed])

            m_data_dict[seed] = m_data

            m_h = plt.scatter(m_x, m_y, color='k', s=7)


            d_h = plt.scatter(d_x, d_y, color='r', s=7, marker='+')
            d_h_forced = plt.scatter(d_x_forced, d_y_forced, color='g', s=7, marker='x')
            d_h_list.append(d_h)
            d_label_list.append('Dale$\chi$; $\chi^2<=\chi^2_{min}+21.03$')
            d_h_list.append(d_h_forced)
            d_label_list.append('Dale$\chi$; $\chi^2<=116.03$')

            text = ('MultiNest: %.2e $\chi^2$ calls\n' % n_mult
                    + 'Dale$\chi$: %.2e $\chi^2$ calls; $\chi^2_{min}=%.2f$' % (limit, chisq_min))
            text_list.append(text)

            for xx in (m_x.min(), m_x.max(), d_x.min(), d_x.max()):
                if xmax is None or xx>xmax:
                    xmax=xx
                if xmin is None or xx<xmin:
                    xmin=xx

            for yy in (m_y.min(), m_y.max(), d_y.min(), d_y.max()):
                if ymax is None or yy>ymax:
                    ymax=yy
                if ymin is None or yy<ymin:
                    ymin=yy
            print "one seed took ",time.time()-t_start

        for i_seed in range(len(seed_list)):
            plt.subplot(3,2,i_seed+1)
            dx=xmax-xmin
            dy=ymax-ymin
            plt.xlim((xmin-0.05*dx, xmax+0.05*dx))
            plt.ylim((ymin-0.05*dy, ymax+0.6*dy))
            if i_seed==0:
                plt.xlabel('$\\theta_%d$' % dim[0], fontsize=15)
                plt.ylabel('$\\theta_%d$' % dim[1], fontsize=15)
            plt.text(xmin, ymax+0.01*dy, text_list[i_seed], fontsize=10)

        plt.legend([m_h] + d_h_list,
                   ['MultiNest'] + d_label_list,
                   fontsize=10,
                   bbox_to_anchor=(1.05,1),
                   loc=2)

        fig_name = 'compare_d12_%d_%d.png' % (dim[0], dim[1])
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, fig_name))
        plt.close()
