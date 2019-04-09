import numpy as np
import os
import time
from scipy import spatial as scipy_spatial

def read_dalex_pts(data_file, raw_dim):
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
            data_pts[i_line] = pt
            chisq[i_line] = params[dim]
            i_line += 1

    assert i_line == n_lines
    return data_pts, chisq

if __name__ == "__main__":

    dalex_file = 'planck_out_high_q2.txt'
    cutoff = 2980.0
    dim=33
    multinest_file = os.path.join('multinest', 'planck_d33_s123_n1000_t1.00e-03.txt')

    dalex_pts, dalex_chisq = read_dalex_pts(dalex_file, dim)
    valid = np.where(dalex_chisq<cutoff)
    dalex_pts = dalex_pts[valid]
    dalex_chisq = dalex_chisq[valid]
    scale = np.zeros(dim,dtype=float)
    target_chisq=dalex_chisq.min()+47.41
    inside = np.where(dalex_chisq<target_chisq)
    inside_pts = dalex_pts[inside]
    for i_dim in range(dim):
        scale[i_dim] = inside_pts[:,i_dim].max()-inside_pts[:,i_dim].min()
    dalex_pts = np.array([vv/scale for vv in dalex_pts])
    print('read dalex data')
    #tree = scipy_spatial.KDTree(dalex_pts, leafsize=10)
    #print('built tree')

    n_multinest = 0
    with open(multinest_file, 'r') as in_file:
        for line in in_file:
            n_multinest += 1

    multinest_file = os.path.join('multinest', 'planck_d33_s123_n1000_t1.00e-03.txt')
    multinest_pts = np.zeros((n_multinest, dim), dtype=float)
    multinest_chisq = np.zeros(n_multinest, dtype=float)
    with open(multinest_file, 'r') as in_file:
        for i_line, line in enumerate(in_file):
            params = np.array(line.strip().split()).astype(float)
            multinest_pts[i_line] = params[2:]/scale
            multinest_chisq[i_line]= params[1]

    to_search = len(multinest_pts)
    print('read multinest data')
    #nn_dist, nn_dex = tree.query(multinest_pts[:to_search], k=1)
    #print('queried')
    t_start = time.time()
    with open('mult_v_dalex.txt', 'w') as out_file:
        out_file.write('# mult dalex dd\n')
        for i_pt in range(to_search):
            if i_pt>0 and i_pt%500 == 0:
                duration = (time.time()-t_start)/3600.0
                predicted = to_search*duration/i_pt
                print('%d of %d -- %.2e, %.2e (hrs)' %
                (i_pt, to_search,duration,predicted))
            dd = np.sqrt(np.sum((multinest_pts[i_pt]-dalex_pts)**2,axis=1))
            min_dex = np.argmin(dd)
            out_file.write('%e %e %e %e\n' %
            (multinest_chisq[i_pt], dalex_chisq[min_dex],
             dd[min_dex],
             multinest_chisq[i_pt]-dalex_chisq[min_dex]))
