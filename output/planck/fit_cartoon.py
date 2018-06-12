import numpy as np
import time
import gc

if __name__ == "__main__":

    order = 7
    dim = 33
    center = np.zeros(dim, dtype=float)
    bases = np.zeros((dim,dim), dtype=float)
    radii = np.zeros(dim, dtype=float)
    basis_file = 'ellipse_bases.txt'
    data_file = 'planck_out_high_q2_culled_cartoon.txt'
    with open(basis_file, 'r') as in_file:
        center_line = in_file.readline()
        params = center_line.strip().split()
        for ii in range(dim):
            mu = float(params[ii])
            center[ii] = mu
        for i_basis, line in enumerate(in_file):
            params = line.strip().split()
            for ii in range(dim):
                mu = float(params[ii])
                bases[i_basis][ii] = mu
            radii[i_basis] = float(params[-1])
            assert len(params) == dim+1

    with open(data_file, 'r') as in_file:
        data_lines = in_file.readlines()

    #data_lines = data_lines[:100000]

    chisq_min = None
    for line in data_lines:
        if line[0]== '#':
            continue
        params = line.strip().split()
        if chisq_min is None or float(params[dim])<chisq_min:
            chisq_min = float(params[dim])
            center = np.array(params[:dim]).astype(float)

    print('got chisq_min %e' % chisq_min)

    n_data = len(data_lines)-1
    n_training = 2*n_data//3

    n_training = 10000

    n_validation = n_data-n_training

    ratio = 1.0-n_training/n_data

    training_pts = np.zeros((n_training, dim), dtype=float)
    training_chisq = np.zeros(n_training, dtype=float)
    training_r = np.zeros(n_training, dtype=float)
    validation_pts = np.zeros((n_validation, dim), dtype=float)
    validation_chisq = np.zeros(n_validation, dtype=float)
    validation_r = np.zeros(n_validation, dtype=float)

    rng = np.random.RandomState(553321984)
    roll = rng.random_sample(n_data)

    set_training = 0
    set_validation = 0
    projected = np.zeros(dim, dtype=float)

    t_start = time.time()
    for i_line, line in enumerate(data_lines[1:]):
        if i_line>0 and i_line%10000==0:
            duration = time.time()-t_start
            predicted = len(data_lines)*duration/i_line
            print(i_line,duration,predicted/60.0,' mins')
        if line[0] == '#':
            continue
        params = np.array(line.strip().split()).astype(float)
        pt = params[:dim]
        for ii in range(dim):
            projected[ii] = np.dot(bases[ii], pt-center)

        rr = np.power(projected/radii,2).sum()

        if (set_validation==n_validation or
            roll[i_line]>ratio and set_training<n_training):

            for ii in range(dim):
                training_pts[set_training][ii] = projected[ii]
            training_chisq[set_training] = params[dim]
            training_r[set_training] = rr
            set_training+=1
        else:
             for ii in range(dim):
                 validation_pts[set_validation][ii] = projected[ii]
             validation_chisq[set_validation] = params[dim]
             validation_r[set_validation] = rr
             set_validation += 1

    del data_lines
    gc.collect()

    assert set_training == n_training
    assert set_validation == n_validation

    n_matrix = order*dim

    #i_matrix = i_dim*order + i_order

    n_iteration = 300

    delta_chisq=47.41
    chisq_min = min(training_chisq.min(), validation_chisq.min())
    sigma_sq = np.ones(len(training_chisq), dtype=float)

    pt_powers = np.zeros(n_training*dim*order, dtype=float)
    # i_power = i_pt*dim*order+i_dim*order+i_order
    # i_power = i_pt*dim*order+i_matrix
    print('made pt_powers')
    for i_pt in range(n_training):
        for i_dim in range(dim):
            for i_order in range(order):
                i_power = i_pt*dim*order+i_dim*order+i_order
                pt_powers[i_power] = training_pts[i_pt][i_dim]**(i_order+1)

    min_failure = None
    for iteration in range(n_iteration):

        mm = np.zeros((n_matrix, n_matrix), dtype=float)
        bb = np.zeros(n_matrix, dtype=float)

        print('starting matrix loop %e' % training_chisq.min())
        t_start = time.time()

        matrix_ct = 0
        for i_matrix_1 in range(n_matrix):
            zz1 = pt_powers[i_matrix_1:len(pt_powers):dim*order]
            mm[i_matrix_1][i_matrix_1] = (zz1**2/sigma_sq).sum()
            bb[i_matrix_1] = ((training_chisq-chisq_min)*zz1/sigma_sq).sum()
            for i_matrix_2 in range(i_matrix_1+1, n_matrix):
                zz2 = pt_powers[i_matrix_2:len(pt_powers):dim*order]
                mu = (zz1*zz2/sigma_sq).sum()
                mm[i_matrix_1][i_matrix_2] = mu
                mm[i_matrix_2][i_matrix_1] = mu
                matrix_ct += 1
                if (matrix_ct+1)%10000 == 0:
                    duration = time.time()-t_start
                    predicted = (0.5*n_matrix*(n_matrix-1))*duration/(matrix_ct+1)
                    print(matrix_ct,duration,predicted/3600.0)

        coeffs = np.linalg.solve(mm, bb)
        assert len(coeffs) == n_matrix

        fit_chisq = np.zeros(n_training, dtype=float)
        for i_pt, pt in enumerate(training_pts):
            vv = 0.0
            for i_dim in range(dim):
                for i_order in range(order):
                    i_matrix = i_dim*order+i_order
                    mu = coeffs[i_matrix]*pt[i_dim]**(i_order+1)
                    vv += mu
            fit_chisq[i_pt] = vv

        chi_wrong = np.abs(training_chisq-chisq_min-fit_chisq)

        # fine
        quad_1 = np.where(np.logical_and(training_chisq-chisq_min>delta_chisq,
                                         fit_chisq>delta_chisq))

        quad_1_1 = np.where(np.logical_and(training_chisq-chisq_min>1.5*delta_chisq,
                                           fit_chisq>1.5*delta_chisq))

        # bad
        quad_2 = np.where(np.logical_and(training_chisq-chisq_min<delta_chisq,
                                         fit_chisq>delta_chisq))

        # fine
        quad_3 = np.where(np.logical_and(training_chisq-chisq_min<delta_chisq,
                                         fit_chisq<delta_chisq))

        # bad
        quad_4 = np.where(np.logical_and(training_chisq-chisq_min>delta_chisq,
                                         fit_chisq<delta_chisq))

        n_offenders = len(quad_4[0])+len(quad_2[0])
        n_non_offenders = len(quad_1[0])+len(quad_3[0])
        max_mismatch = max(chi_wrong[quad_4].max(), chi_wrong[quad_2].max())

        sigma_sq[quad_3] += 2.0

        sigma_sq[quad_1] += 2.0
        sigma_sq[quad_1_1] += 2.0

        max_val = max(chi_wrong[quad_2].max(), chi_wrong[quad_4].max())
        #max_val = chi_wrong[quad_2].max()
        sigma_sq[quad_2] += 1.0-chi_wrong[quad_2]/max_val

        #max_val = chi_wrong[quad_4].max()
        sigma_sq[quad_4] += 1.0-chi_wrong[quad_4]/max_val

        assert sigma_sq.min()>0.99

        if min_failure is None or max_mismatch<min_failure:
            min_failure = max_mismatch
            best_coeffs = coeffs

            with open('coeffs_test.txt', 'w') as out_file:
                for cc in best_coeffs:
                    out_file.write('%e\n' % cc)
            with open('fit_chisq.txt', 'w') as out_file:
                for i_pt in range(n_training):
                    out_file.write('%e %e\n' %
                                   (training_chisq[i_pt]-chisq_min,
                                    fit_chisq[i_pt]))

        print('\niter %d max_mis %.4e best %.4e -- %d %d -- %.2e %.2e %.2e' %
              (iteration, max_mismatch, min_failure, n_offenders, n_non_offenders,
               sigma_sq.min(),np.median(sigma_sq),sigma_sq.max()))

        print('')

    for i_dim in range(dim):
        i_matrix = i_dim*order+order-1
        print(coeffs[i_matrix])
