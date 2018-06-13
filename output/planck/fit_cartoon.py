import numpy as np
import time
import gc


class CartoonFitter(object):

    def __init__(self, order, dim):
        self.order = order
        self.dim = dim

    def read_bases(self, basis_file):
        self.center = np.zeros(self.dim, dtype=float)
        self.bases = np.zeros((self.dim,self.dim), dtype=float)
        self.radii = np.zeros(self.dim, dtype=float)
        with open(basis_file, 'r') as in_file:
            center_line = in_file.readline()
            params = center_line.strip().split()
            for ii in range(self.dim):
                mu = float(params[ii])
                self.center[ii] = mu
            for i_basis, line in enumerate(in_file):
                params = line.strip().split()
                for ii in range(self.dim):
                    mu = float(params[ii])
                    self.bases[i_basis][ii] = mu
                self.radii[i_basis] = float(params[-1])
                assert len(params) == self.dim+1

    def read_data(self, data_file, n_training):

        with open(data_file, 'r') as in_file:
            data_lines = in_file.readlines()

        #data_lines = data_lines[:1000000]

        self.chisq_min = None
        for line in data_lines:
            if line[0]== '#':
                continue
            params = line.strip().split()
            if self.chisq_min is None or float(params[self.dim])<self.chisq_min:
                self.chisq_min = float(params[self.dim])
                self.center = np.array(params[:self.dim]).astype(float)

        print('got chisq_min %e' % self.chisq_min)

        n_data = len(data_lines)-1
        if n_training>1:
            self.n_training = n_training
        else:
            self.n_training = int(np.round(n_training*data_lines))

        self.n_validation = n_data-n_training

        ratio = 1.0-self.n_training/n_data

        self.training_pts = np.zeros((self.n_training, self.dim), dtype=float)
        self.training_chisq = np.zeros(self.n_training, dtype=float)
        self.training_r = np.zeros(self.n_training, dtype=float)
        self.validation_pts = np.zeros((self.n_validation, self.dim), dtype=float)
        self.validation_chisq = np.zeros(self.n_validation, dtype=float)
        self.validation_r = np.zeros(self.n_validation, dtype=float)

        rng = np.random.RandomState(553321984)
        roll = rng.random_sample(n_data)

        set_training = 0
        set_validation = 0
        projected = np.zeros(self.dim, dtype=float)

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
                projected[ii] = np.dot(self.bases[ii], pt-self.center)

            rr = np.power(projected/self.radii,2).sum()

            if (set_validation==self.n_validation or
                roll[i_line]>ratio and set_training<self.n_training):
                self.training_pts[set_training] = projected
                self.training_chisq[set_training] = params[self.dim]
                self.training_r[set_training] = rr
                set_training+=1
            else:
                 self.validation_pts[set_validation] = projected
                 self.validation_chisq[set_validation] = params[self.dim]
                 self.validation_r[set_validation] = rr
                 set_validation += 1

        del data_lines
        gc.collect()

        assert set_training == self.n_training
        assert set_validation == self.n_validation

        self.pt_powers = np.zeros(self.n_training*self.dim*self.order, dtype=float)
        # i_power = i_pt*dim*order+i_dim*order+i_order
        # i_power = i_pt*dim*order+i_matrix
        print('made pt_powers')
        for i_pt in range(self.n_training):
            for i_dim in range(self.dim):
                for i_order in range(self.order):
                    i_power = i_pt*self.dim*self.order+i_dim*self.order+i_order
                    self.pt_powers[i_power] = self.training_pts[i_pt][i_dim]**(i_order+1)

        self.log_training_r = np.log(self.training_r+1.0)
        self.delta_chisq=47.41

    def fit(self, sigma_sq):

        n_matrix = self.order*self.dim+1
        #i_matrix = i_dim*order + i_order

        mm = np.zeros((n_matrix, n_matrix), dtype=float)
        bb = np.zeros(n_matrix, dtype=float)

        print('starting matrix loop')
        t_start = time.time()

        matrix_ct = 0
        for i_matrix_1 in range(n_matrix-1):
            zz1 = self.pt_powers[i_matrix_1:len(self.pt_powers):self.dim*self.order]
            mm[i_matrix_1][i_matrix_1] = (zz1**2/sigma_sq).sum()
            bb[i_matrix_1] = ((self.training_chisq-self.chisq_min)*zz1/sigma_sq).sum()
            for i_matrix_2 in range(i_matrix_1+1, n_matrix-1):
                zz2 = self.pt_powers[i_matrix_2:len(self.pt_powers):self.dim*self.order]
                mu = (zz1*zz2/sigma_sq).sum()
                mm[i_matrix_1][i_matrix_2] = mu
                mm[i_matrix_2][i_matrix_1] = mu
                matrix_ct += 1
                if (matrix_ct+1)%10000 == 0:
                    duration = time.time()-t_start
                    predicted = (0.5*n_matrix*(n_matrix-1))*duration/(matrix_ct+1)
                    print(matrix_ct,duration,predicted/3600.0)

        zz1 = self.log_training_r
        bb[n_matrix-1] = ((self.training_chisq-self.chisq_min)*zz1/sigma_sq).sum()
        mm[n_matrix-1][n_matrix-1] = (zz1**2/sigma_sq).sum()
        for i_matrix_2 in range(n_matrix-1):
            zz2 = self.pt_powers[i_matrix_2:len(self.pt_powers):self.dim*self.order]
            mu = (zz1*zz2/sigma_sq).sum()
            mm[n_matrix-1][i_matrix_2] = mu
            mm[i_matrix_2][n_matrix-1] = mu

        coeffs = np.linalg.solve(mm, bb)
        assert len(coeffs) == n_matrix

        fit_chisq = np.zeros(self.n_training, dtype=float)
        for i_matrix in range(n_matrix-1):
            fit_chisq += coeffs[i_matrix]*self.pt_powers[i_matrix:len(self.pt_powers):self.dim*self.order]
        fit_chisq+=coeffs[n_matrix-1]*self.log_training_r

        chi_wrong = np.abs(self.training_chisq-self.chisq_min-fit_chisq)

        # fine
        quad_1 = np.where(np.logical_and(self.training_chisq-self.chisq_min>self.delta_chisq,
                                         fit_chisq>self.delta_chisq))

        quad_1_1 = np.where(np.logical_and(self.training_chisq-self.chisq_min>1.5*self.delta_chisq,
                                           fit_chisq>1.5*self.delta_chisq))

        # bad
        quad_2 = np.where(np.logical_and(self.training_chisq-self.chisq_min<self.delta_chisq,
                                         fit_chisq>self.delta_chisq))

        # fine
        quad_3 = np.where(np.logical_and(self.training_chisq-self.chisq_min<self.delta_chisq,
                                         fit_chisq<self.delta_chisq))

        quad_3_1 = np.where(np.logical_and(self.training_chisq-self.chisq_min<0.5*self.delta_chisq,
                                           0.5*fit_chisq<self.delta_chisq))


        # bad
        quad_4 = np.where(np.logical_and(self.training_chisq-self.chisq_min>self.delta_chisq,
                                         fit_chisq<self.delta_chisq))

        quad_4_1 = np.where(np.logical_and(self.training_chisq-self.chisq_min>self.delta_chisq,
                            np.logical_and(fit_chisq<self.delta_chisq,
                                           chi_wrong<0.5*self.delta_chisq)))


        output = {}
        output['coeffs'] = coeffs
        output['chi_wrong'] = chi_wrong
        output['quad_1'] = quad_1
        output['quad_2'] = quad_2
        output['quad_3'] = quad_3
        output['quad_4'] = quad_4
        output['quad_1_1'] = quad_1_1
        output['quad_3_1'] = quad_3_1
        output['quad_4_1'] = quad_4_1
        output['fit_chisq'] = fit_chisq
        return output


if __name__ == "__main__":

    order = 7
    dim = 33
    fitter = CartoonFitter(order,dim)
    basis_file = 'ellipse_bases.txt'
    fitter.read_bases(basis_file)
    data_file = 'planck_out_high_q2_culled_cartoon.txt'
    fitter.read_data(data_file, 10000)

    n_iteration = 300
    sigma_sq = np.ones(len(fitter.training_chisq), dtype=float)

    n_matrix = order*dim+1
    #i_matrix = i_dim*order + i_order

    min_failure = None
    for iteration in range(n_iteration):

        results = fitter.fit(sigma_sq)

        chi_wrong = results['chi_wrong']
        coeffs = results['coeffs']
        quad_1 = results['quad_1']
        quad_2 = results['quad_2']
        quad_3 = results['quad_3']
        quad_4 = results['quad_4']
        quad_3_1 = results['quad_3_1']
        quad_1_1 = results['quad_1_1']
        quad_4_meh = results['quad_4_1']
        fit_chisq = results['fit_chisq']

        n_offenders = len(quad_4[0])+len(quad_2[0])
        n_non_offenders = len(quad_1[0])+len(quad_3[0])
        max_mismatch = max(chi_wrong[quad_4].max(), chi_wrong[quad_2].max())

        if min_failure is None or max_mismatch<min_failure:
            min_failure = max_mismatch
            best_coeffs = coeffs

            with open('coeffs_test.txt', 'w') as out_file:
                for cc in best_coeffs:
                    out_file.write('%e\n' % cc)
            with open('fit_chisq.txt', 'w') as out_file:
                for i_pt in range(fitter.n_training):
                    out_file.write('%e %e %e %e\n' %
                                   (fitter.training_chisq[i_pt]-fitter.chisq_min,
                                    fit_chisq[i_pt],
                                    fitter.training_r[i_pt],
                                    fitter.training_chisq[i_pt]-fitter.chisq_min-fit_chisq[i_pt]))

        print('\niter %d max_mis %.4e best %.4e -- %d %d -- %.2e %.2e %.2e' %
              (iteration, max_mismatch, min_failure, n_offenders, n_non_offenders,
               sigma_sq.min(),np.median(sigma_sq),sigma_sq.max()))

        print('')


        sigma_sq[quad_3] += 2.0
        sigma_sq[quad_3_1] += 1.0

        sigma_sq[quad_1] += 2.0
        sigma_sq[quad_1_1] += 2.0

        max_val = max(chi_wrong[quad_2].max(), chi_wrong[quad_4].max())
        #max_val = chi_wrong[quad_2].max()
        sigma_sq[quad_2] += 1.0-chi_wrong[quad_2]/max_val

        #max_val = chi_wrong[quad_4].max()
        sigma_sq[quad_4_meh] += 1.0-chi_wrong[quad_4_meh]/max_val

        assert sigma_sq.min()>0.99

    for i_dim in range(dim):
        i_matrix = i_dim*order+order-1
        print(coeffs[i_matrix])
