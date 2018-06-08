import numpy as np
import time
import gc

if __name__ == "__main__":

    order = 7
    dim = 33
    center = np.zeros(dim, dtype=float)
    bases = np.zeros((dim,dim), dtype=float)
    basis_file = 'ellipse_bases.txt'
    data_file = 'planck_out_high_q2_culled.txt'
    with open(basis_file, 'r') as in_file:
        center_line = in_file.readline()
        params = center_line.strip().split()
        for ii in range(len(params)):
            mu = float(params[ii])
            center[ii] = mu
        for i_basis, line in enumerate(in_file):
            params = line.strip().split()
            for ii in range(len(params)):
                mu = float(params[ii])
                bases[i_basis][ii] = mu

    with open(data_file, 'r') as in_file:
        data_lines = in_file.readlines()
    n_data = len(data_lines)-1
    n_training = 2*n_data//3

    n_training = 3000

    n_validation = n_data-n_training

    ratio = 1.0-n_training/n_data

    training_pts = np.zeros((n_training, dim), dtype=float)
    training_chisq = np.zeros(n_training, dtype=float)
    validation_pts = np.zeros((n_validation, dim), dtype=float)
    validation_chisq = np.zeros(n_validation, dtype=float)

    rng = np.random.RandomState(553321984)
    roll = rng.random_sample(n_data)

    set_training = 0
    set_validation = 0
    projected = np.zeros(dim, dtype=float)
    for i_line, line in enumerate(data_lines[1:]):
        if line[0] == '#':
            continue
        params = np.array(line.strip().split()).astype(float)
        pt = params[:dim]
        for ii in range(dim):
            projected[ii] = np.dot(bases[ii], pt-center)

        if (set_validation==n_validation or
            roll[i_line]>ratio and set_training<n_training):

            for ii in range(dim):
                training_pts[set_training][ii] = projected[ii]
            training_chisq[set_training] = params[dim]
            set_training+=1
        else:
             for ii in range(dim):
                 validation_pts[set_validation][ii] = projected[ii]
             validation_chisq[set_validation] = params[dim]
             set_validation += 1

    del data_lines
    gc.collect()

    assert set_training == n_training
    assert set_validation == n_validation

    n_matrix = order*dim

    mm = np.zeros((n_matrix, n_matrix), dtype=float)
    bb = np.zeros(n_matrix, dtype=float)

    #i_matrix = i_dim*order + i_order

    delta_chisq=47.41
    chisq_min = min(training_chisq.min(), validation_chisq.min())
    sigma_sq = np.where(training_chisq<chisq_min+delta_chisq,
                        1.0,
                        ((training_chisq-chisq_min)/delta_chisq)**2)

    print('starting matrix loop')
    t_start = time.time()
    for i_pt in range(n_training):
        cc = training_chisq[i_pt]-chisq_min
        ss = sigma_sq[i_pt]
        for i_matrix_1 in range(n_matrix):
            i_dim_1 = i_matrix_1//order
            i_order_1 = i_matrix_1%order
            zz1 = training_pts[i_pt][i_dim_1]**(i_order_1+1)
            mm[i_matrix_1][i_matrix_1] += zz1**2/ss
            bb[i_matrix_1] += cc*zz1/ss
            for i_matrix_2 in range(i_matrix_1+1, n_matrix):
                i_dim_2 = i_matrix_2//order
                i_order_2 = i_matrix_2%order
                zz2 = training_pts[i_pt][i_dim_2]**(i_order_2+1)
                mm[i_matrix_1][i_matrix_2] += zz1*zz2/ss
                mm[i_matrix_2][i_matrix_1] += zz1*zz2/ss

        if (i_pt+1)%1000 == 0:
            duration = time.time()-t_start
            predicted = n_training*duration/(i_pt+1)
            print(i_pt,duration,predicted/3600.0)

    print('solving coeffs')
    coeffs = np.linalg.solve(mm, bb)
    assert len(coeffs) == n_matrix
    print(coeffs.max(),coeffs.min(),np.median(coeffs))

    print('writing junk.txt')
    with open("junk.txt", "w") as out_file:
        for i_pt, pt in enumerate(training_pts):
            vv = 0.0
            for i_dim in range(dim):
                for i_order in range(order):
                    i_matrix = i_dim*order+i_order
                    mu = coeffs[i_matrix]*pt[i_dim]**(i_order+1)
                    vv += mu
            out_file.write('%e %e %e\n' %
                           (training_chisq[i_pt]-chisq_min, vv,
                            training_chisq[i_pt]-chisq_min-vv))
