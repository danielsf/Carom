import numpy as np
import os
import argparse

def find_basis(pts, fn):
    min_dex = np.argmin(fn)
    min_pt = pts[min_dex]
    xx = pts.transpose()
    dim = len(xx)

    covar = np.array([((xx[ii]-min_pt[ii])*(fn[ii]-fn[min_dex])).sum()
                      for ii in range(dim)])
    covar /= float(len(pts))
    return covar


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', type=str, default='.')
    args = parser.parse_args()

    dim = 12
    dt_list = []
    for i in range(dim):
        dt_list.append(('p%d' % i, float))
    dt_list.append(('chisq', float))
    dt_list.append(('dex', int))
    data_dtype = np.dtype(dt_list)
    end_pt_file = os.path.join(args.dir, 'output_delta.txt_end_pts.txt')
    data = np.genfromtxt(end_pt_file, dtype=data_dtype)

    bases_dtype = np.dtype([('e%d' % ii, float) for ii in range(dim)])
    bases_file = os.path.join(args.dir, 'quad_bases.txt')
    bases =np.genfromtxt(bases_file, dtype=bases_dtype).transpose()

    fn = data['chisq']
    min_dex = np.argmin(fn)
    fn_min = fn[min_dex]
    pts = np.array([data['p%d'%ii] for ii in range(dim)]).transpose()
    min_pt = pts[min_dex]
    pts -= min_pt
    pts_proj = np.zeros(pts.shape, dtype=float)
    for ii in range(len(pts)):
        for jj in range(dim):
           pts_proj[ii][jj]=0.0
           for kk in range(dim):
               pts_proj[ii][jj] += bases[jj][kk]*pts[ii][kk]

    xx = pts_proj.transpose()

    sigma = np.zeros(len(pts), dtype=float)
    wgts_file = os.path.join(args.dir, 'quad_wgts.txt')
    with open(wgts_file, 'r') as in_file:
        for i_line, line in enumerate(in_file):
            ss = float(line.strip())
            sigma[i_line]=ss
    sigma_sq = sigma**2

    mm = np.zeros((dim, dim), dtype=float)
    bb = np.zeros(dim, dtype=float)
    for ii in range(dim):
        bb[ii] = ((xx[ii]**2)*(fn-fn_min)/sigma_sq).sum()
        for jj in range(dim):
            mm[ii][jj] = ((xx[ii]**2)*(xx[jj]**2)/sigma_sq).sum()

    aa = np.linalg.solve(mm, bb)

    test_file = os.path.join(args.dir, 'quadratic_test.txt')
    with open(test_file, 'w') as out_file:
        for i_pt in range(len(pts)):
            chisq_fit = fn_min + (aa*(pts_proj[i_pt]**2)).sum()
            out_file.write('%e %e %e - %e %e -- ' %
                           (fn[i_pt], chisq_fit, np.abs(fn[i_pt]-chisq_fit),
                            (pts[i_pt]**2).sum(), pts_proj[i_pt][5]))

            for ii in range(dim):
                out_file.write('%e ' % pts_proj[i_pt][ii])
            out_file.write('\n')

    a_wgts_file = os.path.join(args.dir, 'quad_a_wgts.txt')
    with open(a_wgts_file, 'w') as out_file:
        for aval in aa:
            out_file.write('%e\n' % aval)

    print(aa)
    print(min_pt)
