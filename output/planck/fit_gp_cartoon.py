import os
import numpy as np

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

def get_gp_set(data_pts, chisq, radii, n_pts, dim):
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
            valid = np.where(chisq<target_chisq+tol)
            pts_considered = normalized_pts[valid]
            dd_min = None
        elif ii == n_pts//2:
            valid = np.where(np.abs(chisq-target_chisq-4.0*tol)<tol)
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
        if dd_min is None:
            dd_min = dd
        else:
            dd_min = np.where(dd_min<dd, dd_min, dd)

        chosen_pt = np.argmax(dd_min)
        actual_dex = valid[0][chosen_pt]
        #assert actual_dex not in already_chosen
        #already_chosen.add(actual_dex)
        gp_pts[n_assigned] = data_pts[actual_dex]
        gp_chisq[n_assigned] = chisq[actual_dex]
        n_assigned += 1
        print(n_assigned,dd_min[chosen_pt])

    return gp_pts, gp_chisq

if __name__ == "__main__":

    dim = 33
    basis_file = 'ellipse_bases.txt'
    data_file = 'planck_out_high_q2_culled_cartoon.txt'

    bases, radii = read_bases(basis_file, dim)
    data_pts, chisq = project_pts(data_file, bases, dim)

    print(chisq.min(),np.median(chisq),chisq.max())

    gp_pts, gp_chisq = get_gp_set(data_pts, chisq, radii, 2000, dim)

    with open('gp_pts.txt', 'w') as out_file:
        for i_pt in range(len(gp_pts)):
            for ii in range(dim):
                out_file.write('%e ' % gp_pts[i_pt][ii])
            out_file.write('%e\n' % gp_chisq[i_pt])
