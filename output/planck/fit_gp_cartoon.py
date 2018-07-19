import os
import numpy as np
import time
from scipy import spatial as scipy_spatial
from iminuit import Minuit

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make_2d_histogram(xx, yy, dx, dy):
    """
    returns indices and counts of unique points on the map
    """
    i_color1 = np.round(xx/dx).astype(int)
    i_color2 = np.round(yy/dy).astype(int)
    dex_reverse = np.array([i_color1, i_color2])
    dex_arr = dex_reverse.transpose()
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    dex_raw = np.ascontiguousarray(dex_arr).view(np.dtype((np.void, dex_arr.dtype.itemsize*dex_arr.shape[1])))
    _, unique_rows, unique_counts = np.unique(dex_raw, return_index=True, return_counts=True)

    return unique_rows, unique_counts


def plot_color(xx, yy, dx, dy):
    dexes, cts = make_2d_histogram(xx, yy, dx, dy)
    sorted_dex = np.argsort(cts)
    dexes = dexes[sorted_dex]
    cts = cts[sorted_dex]
    plt.scatter(xx[dexes], yy[dexes], c=cts, s=5,
                cmap=plt.cm.gist_ncar, edgecolor='')

    plt.colorbar()


def plot_color_mesh(xx, yy, dx, dy, vmin=None, vmax=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex')

    x_mesh=np.arange(xx.min(),xx.max()+dx,dx)
    y_mesh=np.arange(yy.min(),yy.max()+dy,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    z_mesh = np.zeros(shape=x_mesh.shape, dtype=int)
    ct_1000b = 0

    for dex, ct in zip(dex_list, ct_list):
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]
        z_mesh[iy][ix] += ct

    z_mesh = np.ma.masked_where(z_mesh==0,z_mesh)
    z_mesh = z_mesh/z_mesh.sum()
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    plt.colorbar()


def project_pts(data_file, bases, raw_dim):
    n_lines = 0
    with open(data_file, 'r') as in_file:
        for line in in_file:
            if not line.startswith('#'):
                n_lines += 1

    data_pts = np.zeros((n_lines, dim), dtype=float)
    chisq = np.zeros(n_lines, dtype=float)

    i_line = 0
    with open(data_file, 'r') as in_file:
        for line in in_file:
            if line.startswith('#'):
                continue
            params = np.array(line.strip().split()).astype(float)
            pt = params[:raw_dim]
            projected = np.dot(bases, pt)
            data_pts[i_line] = projected
            chisq[i_line] = params[dim]
            i_line += 1

    assert i_line == n_lines
    return data_pts, chisq

def read_bases(basis_file, raw_dim):
    center = np.zeros(raw_dim, dtype=float)
    bases = np.zeros((raw_dim,raw_dim), dtype=float)
    radii = np.zeros(raw_dim, dtype=float)
    with open(basis_file, 'r') as in_file:
        center_line = in_file.readline()
        params = center_line.strip().split()
        for ii in range(raw_dim):
            mu = float(params[ii])
            center[ii] = mu
        for i_basis, line in enumerate(in_file):
            params = line.strip().split()
            for ii in range(raw_dim):
                mu = float(params[ii])
                bases[i_basis][ii] = mu
            radii[i_basis] = float(params[-1])
    return bases, radii

def get_gp_set(data_pts, chisq, quad_chisq, radii, n_pts, dim):
    delta_chisq = 47.41

    normalized_pts = np.array([vv/radii for vv in data_pts])

    min_dex = np.argmin(chisq)
    center = data_pts[min_dex]
    chisq_min = chisq[min_dex]
    target_chisq = chisq_min+delta_chisq
    print('target %e' % target_chisq)

    gp_pts = np.zeros((n_pts, dim), dtype=float)
    gp_chisq = np.zeros(n_pts, dtype=float)

    gp_pts[0] = center
    gp_chisq[0] = chisq[min_dex]

    tol=0.1*delta_chisq
    n_assigned = 1
    already_chosen = set()

    for ii in range(n_pts-1):
        if ii == 0:
            valid = np.where(np.logical_and(chisq<target_chisq,
                                            quad_chisq<target_chisq))
            pts_considered = normalized_pts[valid]
            dd_min = None
        elif ii==n_pts//10:
            valid = np.where(np.logical_and(quad_chisq>target_chisq, chisq<target_chisq))
            pts_considered = normalized_pts[valid]
            dd_min = None
            for i_old in range(n_assigned):
                dd = np.sum((pts_considered-gp_pts[i_old]/radii)**2, axis=1)
                if dd_min is None:
                    dd_min = dd
                else:
                    dd_min = np.where(dd_min<dd, dd_min, dd)
        elif ii == 3*n_pts//4:
            valid = np.where(np.logical_and(np.abs(chisq-target_chisq-4.0*tol)<tol,
                                            quad_chisq>target_chisq))
            pts_considered = normalized_pts[valid]
            dd_min = None
            for i_old in range(n_assigned):
                dd = np.sum((pts_considered-gp_pts[i_old]/radii)**2, axis=1)
                if dd_min is None:
                    dd_min = dd
                else:
                    dd_min = np.where(dd_min<dd, dd_min, dd)

        dd = np.sum((pts_considered-gp_pts[n_assigned-1]/radii)**2, axis=1)
        #assert len(dd) == len(pts_considered)
        n_update=0
        if dd_min is None:
            dd_min = dd
        else:
            update = np.where(dd<dd_min)
            n_update = len(update[0])
            dd_min[update] = dd[update]

        chosen_pt = np.argmax(dd_min)
        actual_dex = valid[0][chosen_pt]
        #assert actual_dex not in already_chosen
        #already_chosen.add(actual_dex)
        gp_pts[n_assigned] = data_pts[actual_dex]
        gp_chisq[n_assigned] = chisq[actual_dex]
        n_assigned += 1
        print(n_assigned,dd_min[chosen_pt],n_update)

    return gp_pts, gp_chisq


class MultinestMinimizer(object):

    def __init__(self, basis_file=None, data_file=None, multinest_file=None):
        self.dim = 33
        self._metric_min=3.0e30

        self.bases, self._baseline_radii = read_bases(basis_file, self.dim)
        self.data_pts, self.chisq = project_pts(data_file, self.bases, self.dim)


        self.min_dex = np.argmin(self.chisq)
        self.min_pt = np.copy(self.data_pts[self.min_dex])
        self.chisq_min = self.chisq[self.min_dex]

        for ipt in range(len(self.data_pts)):
            self.data_pts[ipt] -= self.min_pt
        self.data_pts_sq = self.data_pts**2

        print(self.chisq.min(),np.median(self.chisq),self.chisq.max())

        ### read in multinest
        n_multinest = 0
        with open(multinest_file, 'r') as in_file:
            for line in in_file:
                n_multinest += 1

        self.multinest_pts = np.zeros((n_multinest, dim), dtype=float)
        self.multinest_chisq = np.zeros(n_multinest, dtype=float)
        with open(multinest_file, 'r') as in_file:
            for i_line, line in enumerate(in_file):
                params = np.array(line.strip().split()).astype(float)
                pt = np.dot(self.bases, params[2:])
                self.multinest_pts[i_line] = pt - self.min_pt
                self.multinest_chisq[i_line]= params[1]

        self.multinest_pts_sq = self.multinest_pts**2
        self._baseline_metric = self.fit_mean_model(self._baseline_radii)

    def _set_rr(self, radii):

        wgts = (1.0/radii)**2
        rr_arr = np.sqrt(np.dot(self.data_pts_sq, wgts))
        multinest_rr = np.sqrt(np.dot(self.multinest_pts_sq, wgts))

        return rr_arr, multinest_rr

    def fit_mean_model(self, radii):
        t_start = time.time()

        (rr_arr,
         multinest_rr) = self._set_rr(radii)

        #print('set_rr took %e seconds' % (time.time()-t_start))

        sorted_dex = np.argsort(rr_arr)
        rr_arr = rr_arr[sorted_dex]
        chisq_sorted = self.chisq[sorted_dex]

        #### fit mean model
        n_steps = 4
        d_n = len(rr_arr)//n_steps
        chisq_rr_grid = [self.chisq_min]
        rr_grid = [0.0]
        rr_min = 0.0
        for i_step in range(n_steps):
            i_min = i_step*d_n
            i_max = i_min + d_n
            if i_step == n_steps-1:
                i_max = -1
            rr = 0.5*(rr_arr[i_min]+rr_arr[i_max])
            cc = np.median(chisq_sorted[i_min:i_max])
            rr_grid.append(rr)
            chisq_rr_grid.append(cc)

        rr_grid = np.array(rr_grid)
        chisq_rr_grid = np.array(chisq_rr_grid)
        sorted_dex = np.argsort(rr_grid)
        rr_grid = rr_grid[sorted_dex]
        chisq_rr_grid = chisq_rr_grid[sorted_dex]

        mean_chisq = np.interp(multinest_rr, rr_grid, chisq_rr_grid)

        target = self.chisq_min+47.41
        mis_char = np.where(np.logical_or(
                            np.logical_and(mean_chisq<target, self.multinest_chisq>target),
                            np.logical_and(mean_chisq>target, self.multinest_chisq<target)))

        wgts = np.zeros(len(mean_chisq), dtype=float)
        wgts[mis_char] = 1.0

        metric_terms = np.abs(mean_chisq-self.multinest_chisq)
        max_term = metric_terms[mis_char].max()
        med_term = np.median(metric_terms[mis_char])

        metric = np.sum((metric_terms*wgts)**2)
        if metric<self._metric_min:
            self._metric_min=metric
        if hasattr(self, '_baseline_metric'):
            print('metric %.4e -- %.4e -- %.4e %d %.4e %.4e' %
            (metric,self._metric_min,self._baseline_metric, len(mis_char[0]),
             max_term,med_term))
        else:
            print("metric %e" % metric)
        #print('metric took %e seconds %e' % (time.time()-t_start,metric))
        return metric


if __name__ == "__main__":

    dim = 33
    delta_chisq=47.41
    basis_file = 'ellipse_bases.txt'
    data_file = 'planck_out_high_q2.txt'
    multinest_file = os.path.join('multinest', 'planck_d33_s123_n1000_t1.00e-03.txt')

    rng = np.random.RandomState(57623)

    fitter = MultinestMinimizer(data_file=data_file,
                                multinest_file=multinest_file,
                                basis_file=basis_file)

    valid = np.where(fitter.chisq<fitter.chisq.min()+delta_chisq)
    print('valid %d' % len(valid[0]))

    pp = []
    init_dict = {}
    for ii in range(33):
        pp.append('p%d' % ii)
        init_dict['p%d'%ii] = 0.5*(fitter.data_pts[valid,ii].max()-fitter.data_pts[valid,ii].min())
        init_dict['error_p%d'%ii] = 0.1*init_dict['p%d'%ii]

    mm = Minuit(fitter.fit_mean_model, use_array_call=True, forced_parameters=pp ,**init_dict)
    mm.migrad()
    exit()

    plt.figure(figsize=(30,30))
    valid = np.where(np.logical_and(multinest_chisq<5000.0,mean_chisq<5000.0))
    print('valid points %.5e of %.5e' % (len(valid[0]),len(mean_chisq)))
    plot_color_mesh(multinest_chisq[valid],mean_chisq[valid],1.0,1.0)
    plt.xlabel('chisq', fontsize=40)
    plt.ylabel('mean model', fontsize=40)
    c_max = max(multinest_chisq[valid].max(),mean_chisq[valid].max())
    c_min = min(multinest_chisq[valid].min(),mean_chisq[valid].min())
    #plt.plot([c_min,c_max],[c_min,c_max],linestyle='--',color='r',linewidth=5)
    plt.xlim((multinest_chisq[valid].min(),multinest_chisq[valid].max()))
    plt.ylim((mean_chisq[valid].min(),mean_chisq[valid].max()))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig('mult_mean_density.png')
    exit()


    #### find ell interp grids
    #target_chisq = chisq.min()+target_chisq
    #valid = np.where(chisq<target_chisq+0.1*delta_chisq)
    #ell_pts = np.array([vv/radii for vv in data_pts[valid]])
    #ell_chisq = chisq[valid]
    #ell_tree = scipy_spatial.KDTree(ell_pts, leafsize=100)

    #ell_grid = np.arange(0.1, 2.0, 0.1)
    #prev_pairs = None

    ### select points for GP ###
    mean_chisq = np.interp(rr_arr, rr_grid, chisq_rr_grid)
    n_gp_pts = 2000
    gp_pts_file = 'gp_pts.txt'
    if not os.path.exists(gp_pts_file):
        gp_pts, gp_chisq = get_gp_set(data_pts, chisq, mean_chisq, radii, n_gp_pts, dim)

        with open(gp_pts_file, 'w') as out_file:
            for i_pt in range(len(gp_pts)):
                for ii in range(dim):
                    out_file.write('%e ' % gp_pts[i_pt][ii])
                out_file.write('%e\n' % gp_chisq[i_pt])
    else:
        print('reading in gp_pts')
        gp_pts = np.zeros((n_gp_pts, dim), dtype=float)
        gp_chisq = np.zeros(n_gp_pts, dtype=float)
        with open(gp_pts_file, 'r') as in_file:
            for i_line, line in enumerate(in_file):
                params = np.array(line.strip().split()).astype(float)
                gp_pts[i_line] = params[:dim]
                gp_chisq[i_line] = params[dim]
        assert i_line == n_gp_pts-1

    #### build GP covariance matrix
    normalized_pts = np.array([vv/radii for vv in gp_pts])

    tree = scipy_spatial.KDTree(normalized_pts, leafsize=100)

    n_neighbors = 10

    covar_dist, covar_dex = tree.query(normalized_pts, k=n_neighbors)

    print('covar_dex ',covar_dex.shape)

    covar_matrix = -1.0*np.ones((len(gp_pts), len(gp_pts)), dtype=float)

    covar_weights = np.arange(1.0, 0.0, -1.0/n_neighbors)
    assert covar_weights.min()>0.0

    print('building covar')
    for i_1 in range(len(covar_dex)):
        ell = covar_dist[i_1][1]
        ddsq = np.sum(((gp_pts[i_1]-gp_pts)/radii)**2, axis=1)
        local_wgt = np.exp(-0.5*ddsq/(ell*ell))
        for i_2 in range(n_gp_pts):
            ww = local_wgt[i_2]
            if covar_matrix[i_1][i_2]<-0.1 or ww<covar_matrix[i_1][i_2]:
                covar_matrix[i_1][i_2] = ww
                covar_matrix[i_2][i_1] = ww

    untouched = np.where(covar_matrix<-0.1)
    covar_matrix[untouched] = 0.0
    print('built covar %e %e %e' %
    (covar_matrix.min(),np.median(covar_matrix),covar_matrix.max()))

    nugget = 1.0e-5
    for i_1 in range(len(covar_matrix)):
        covar_matrix[i_1][i_1] += nugget

    print('inverting covar')
    covar_inv = np.linalg.inv(covar_matrix)

    rr_gp = np.zeros(n_gp_pts, dtype=float)
    rr_gp = np.sqrt(np.sum(((gp_pts-min_pt)/radii)**2,axis=1))
    gp_quad_chisq = np.interp(rr_gp,rr_grid,chisq_rr_grid)

    covar_vec = np.dot(covar_inv, gp_chisq-gp_quad_chisq)
    assert len(covar_vec) == n_gp_pts

    #### test against MultiNest samples

    with open('multinest_gp_comparison.txt', 'w') as out_file:
        out_file.write('# true gp quad\n')
        for pp, cc, in zip(multinest_pts, multinest_chisq):
            rr = np.sqrt(np.sum(((pp-min_pt)/radii)**2))
            qq = np.interp(rr, rr_grid, chisq_rr_grid)
            normed = pp/radii
            dist, dex = tree.query(normed, k=n_neighbors)
            ell = dist[1]
            ddsq = np.sum(((pp-gp_pts)/radii)**2, axis=1)
            cq = np.exp(-0.5*ddsq/(ell*ell))
            #cq = np.zeros(n_gp_pts, dtype=float)
            #cq[dex] = covar_weights
            fit = qq + np.dot(cq, covar_vec)

            out_file.write('%e %e %e %e\n' % (cc, fit, qq, rr))
    exit()

    with open('multinest_gp_comparison.txt', 'w') as out_file:
        out_file.write('# true gp quad\n')
        for pp, cc, cq_cheat in zip(gp_pts, gp_chisq, covar_matrix):
            rr = np.sqrt(np.sum(((pp-min_pt)/radii)**2))
            qq = np.interp(rr, rr_grid, chisq_rr_grid)
            normed = pp/radii
            dist, dex = tree.query(normed, k=n_neighbors)
            ell = dist[1]
            ddsq = np.sum(((pp-gp_pts)/radii)**2, axis=1)
            cq = np.exp(-0.5*ddsq/(ell*ell))
            #cq = np.zeros(n_gp_pts, dtype=float)
            #cq[dex] = covar_weights
            dd = np.sqrt(np.sum((cq-cq_cheat)**2))
            fit = qq + np.dot(cq, covar_vec)

            out_file.write('%e %e %e -- %e\n' % (cc, fit, qq, dd))
