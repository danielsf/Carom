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
    delta_chisq=47.41
    basis_file = 'ellipse_bases.txt'
    data_file = 'planck_out_high_q2.txt'

    rng = np.random.RandomState(57623)

    bases, radii = read_bases(basis_file, dim)
    data_pts, chisq = project_pts(data_file, bases, dim)

    min_dex = np.argmin(chisq)
    min_pt = data_pts[min_dex]
    chisq_min = chisq[min_dex]

    print(chisq.min(),np.median(chisq),chisq.max())

    n_gp_pts = 2000
    gp_pts_file = 'gp_pts.txt'
    if not os.path.exists(gp_pts_file):
        gp_pts, gp_chisq = get_gp_set(data_pts, chisq, radii, n_gp_pts, dim)

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

    #### fit quadratic mean model
    bb = np.zeros(dim, dtype=float)
    mm = np.zeros((dim,dim), dtype=float)
    valid = np.where(np.logical_and(chisq>chisq_min+delta_chisq,chisq<1.0e5))
    quad_pts = data_pts[valid]
    quad_chisq = chisq[valid]

    d_chisq = quad_chisq-chisq_min
    wgt = np.log10(d_chisq)
    for i_dim in range(dim):
        print('idim %d %d' % (i_dim, len(quad_chisq)))
        z1 = ((quad_pts[:,i_dim]-min_pt[i_dim])/radii[i_dim])**2
        bb[i_dim] = np.sum(wgt*d_chisq*z1)
        for i_dim_2 in range(i_dim, dim):
            z2 = ((quad_pts[:,i_dim_2]-min_pt[i_dim_2])/radii[i_dim])**2
            mm[i_dim][i_dim_2] = np.sum(wgt*z1*z2)
            if i_dim != i_dim_2:
                mm[i_dim_2][i_dim] = mm[i_dim][i_dim_2]

    quadratic_coeffs = np.linalg.solve(mm, bb)
    for qq in quadratic_coeffs:
        print(qq)

    metric = 0.0
    with open('quadratic_check.txt', 'w') as out_file:
        for i_pt in range(len(data_pts)):
            fit = np.sum(quadratic_coeffs*((data_pts[i_pt]-min_pt)/radii)**2)
            true_d = chisq[i_pt]-chisq_min
            metric += (true_d-fit)**2
            if ((fit<delta_chisq and true_d>delta_chisq) or
                (fit>delta_chisq and true_d<delta_chisq)):

                out_file.write('%e %e\n' % (chisq[i_pt]-chisq_min, fit))

    pos_coeffs = quadratic_coeffs[np.where(quadratic_coeffs>0.0)]
    filler_coeff = np.median(pos_coeffs)
    quadratic_coeffs = np.where(quadratic_coeffs>0.0, quadratic_coeffs, filler_coeff)
    pos_metric = 0.0
    with open('quadratic_check_pos.txt', 'w') as out_file:
        for i_pt in range(len(data_pts)):
            fit = np.sum(quadratic_coeffs*((data_pts[i_pt]-min_pt)/radii)**2)
            true_d = chisq[i_pt]-chisq_min
            pos_metric += (true_d-fit)**2
            if ((fit<delta_chisq and true_d>delta_chisq) or
                (fit>delta_chisq and true_d<delta_chisq)):

                out_file.write('%e %e\n' % (chisq[i_pt]-chisq_min, fit))

    print('metric %e pos %e\n' % (metric,pos_metric))
    print('filler %e\n' % filler_coeff)
