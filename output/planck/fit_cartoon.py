import numpy as np
import time
import gc
import multiprocessing as mproc
import argparse

def make_matrix(i_proc, pt_powers, dim, order,
                training_chisq, sigma_sq, out_dict):
    t_start = time.time()
    n_matrix = order*dim
    mm = np.zeros((n_matrix, n_matrix), dtype=float)
    bb = np.zeros(n_matrix, dtype=float)
    matrix_ct = 0
    for i_matrix_1 in range(n_matrix):
        zz1 = pt_powers[i_matrix_1:len(pt_powers):dim*order]
        mm[i_matrix_1][i_matrix_1] = (zz1**2/sigma_sq).sum()
        bb[i_matrix_1] = (training_chisq*zz1/sigma_sq).sum()
        for i_matrix_2 in range(i_matrix_1+1, n_matrix):
            zz2 = pt_powers[i_matrix_2:len(pt_powers):dim*order]
            mu = (zz1*zz2/sigma_sq).sum()
            mm[i_matrix_1][i_matrix_2] = mu
            mm[i_matrix_2][i_matrix_1] = mu
            matrix_ct += 1
            if (matrix_ct+1)%10000 == 0:
                duration = time.time()-t_start
                predicted = (0.5*n_matrix*(n_matrix-1))*duration/(matrix_ct+1)
                print(matrix_ct,duration,' sec ',predicted/60.0,' mins')

    out_dict['m%d' % i_proc]= mm
    out_dict['b%d' % i_proc] = bb
    return None


class CartoonFitter(object):

    def __init__(self, order, dim):
        self.order = order
        self.dim = dim+1
        self.raw_dim = dim
        self.i_iteration = 0

    def read_bases(self, basis_file):
        self.center = np.zeros(self.raw_dim, dtype=float)
        self.bases = np.zeros((self.raw_dim,self.raw_dim), dtype=float)
        self.radii = np.zeros(self.raw_dim, dtype=float)
        with open(basis_file, 'r') as in_file:
            center_line = in_file.readline()
            params = center_line.strip().split()
            for ii in range(self.raw_dim):
                mu = float(params[ii])
                self.center[ii] = mu
            for i_basis, line in enumerate(in_file):
                params = line.strip().split()
                for ii in range(self.raw_dim):
                    mu = float(params[ii])
                    self.bases[i_basis][ii] = mu
                self.radii[i_basis] = float(params[-1])
                assert len(params) == self.raw_dim+1

    def read_data(self, data_file, n_training, n_procs):

        with open(data_file, 'r') as in_file:
            data_lines = in_file.readlines()

        #data_lines = data_lines[:1000000]

        self.chisq_min = None
        for line in data_lines:
            if line[0]== '#':
                continue
            params = line.strip().split()
            if self.chisq_min is None or float(params[self.raw_dim])<self.chisq_min:
                self.chisq_min = float(params[self.raw_dim])
                self.center = np.array(params[:self.raw_dim]).astype(float)

        print('got chisq_min %e' % self.chisq_min)

        n_data = len(data_lines)-1
        if n_training>1:
            self.n_training = int(n_training)
        else:
            self.n_training = int(np.round(n_training*n_data))

        self.n_validation = n_data-self.n_training

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
        projected = np.zeros(self.raw_dim, dtype=float)

        t_start = time.time()
        for i_line, line in enumerate(data_lines[1:]):
            if i_line>0 and i_line%10000==0:
                duration = time.time()-t_start
                predicted = len(data_lines)*duration/i_line
                print(i_line,duration,predicted/60.0,' mins')
            if line[0] == '#':
                continue
            params = np.array(line.strip().split()).astype(float)
            pt = params[:self.raw_dim]
            for ii in range(dim):
                projected[ii] = np.dot(self.bases[ii], pt-self.center)

            rr = np.power(projected/self.radii,2).sum()

            if (set_validation==self.n_validation or
                roll[i_line]>ratio and set_training<self.n_training):
                self.training_pts[set_training][:self.raw_dim] = projected
                self.training_pts[set_training][self.dim-1] = np.log(1.0+rr)
                self.training_chisq[set_training] = params[self.raw_dim]
                self.training_r[set_training] = rr
                set_training+=1
            else:
                 self.validation_pts[set_validation][:self.raw_dim] = projected
                 self.validation_pts[set_validation][self.dim-1] = np.log(1.0+rr)
                 self.validation_chisq[set_validation] = params[self.raw_dim]
                 self.validation_r[set_validation] = rr
                 set_validation += 1

        del data_lines
        gc.collect()

        assert set_training == self.n_training
        assert set_validation == self.n_validation

        self.n_procs = n_procs
        self.pt_powers = {}
        self.pt_dexes = {}
        n_pts = [self.n_training//n_procs]*n_procs
        if sum(n_pts)<self.n_training:
            n_pts[-1] += self.n_training - sum(n_pts)

        # i_power = i_pt*dim*order+i_dim*order+i_order
        # i_power = i_pt*dim*order+i_matrix

        i_start = 0
        pt_sum = 0
        for i_proc in range(n_procs):
            i_end = i_start + n_pts[i_proc]
            self.pt_dexes[i_proc] = (i_start, i_end)
            local_pt_powers = np.zeros((i_end-i_start)*self.dim*self.order, dtype=float)
            for i_pt in range(i_start, i_end):
                pt_sum += 1
                for i_order in range(self.order):
                    #i_power = (i_pt-i_start)*self.dim*self.order+i_dim*self.order+i_order
                    i_power_0 = (i_pt-i_start)*self.dim*self.order+i_order
                    i_power_1 = (i_pt-i_start)*self.dim*self.order+self.dim*order+i_order
                    local_pt_powers[i_power_0:i_power_1:order] = self.training_pts[i_pt]**(i_order+1)
            self.pt_powers[i_proc] = local_pt_powers
            i_start = i_end

        assert pt_sum == self.n_training
        print('made pt_powers')

        self.log_training_r = np.log(self.training_r+1.0)
        self.delta_chisq=47.41
        self.best_mismatch = 2.0e30

    def fit(self, sigma_sq):
        self.i_iteration += 1

        n_matrix = self.order*self.dim
        #i_matrix = i_dim*order + i_order

        #def make_matrix(i_proc, pt_powers, n_matrix, dim, order,
        #        training_chisq, chisq_min, sigma_sq, out_dict):

        print('starting matrix loop')
        mm_dict = {}
        if self.n_procs == 1:
            make_matrix(0, self.pt_powers[0], self.dim, self.order,
                        self.training_chisq-self.chisq_min, sigma_sq,
                        mm_dict)

            mm = mm_dict['m0']
            bb = mm_dict['b0']
        else:
            p_list = []
            mm_dict = mproc.Manager().dict()
            for i_proc in range(self.n_procs):
                i_start = self.pt_dexes[i_proc][0]
                i_end = self.pt_dexes[i_proc][1]
                p = mproc.Process(target=make_matrix,
                                  args=(i_proc, self.pt_powers[i_proc],
                                        self.dim, self.order,
                                        self.training_chisq[i_start:i_end]-self.chisq_min,
                                        sigma_sq[i_start:i_end],
                                        mm_dict))
                p.start()
                p_list.append(p)
            for p in p_list:
                p.join()

            mm = mm_dict['m0']
            bb = mm_dict['b0']
            for i_proc in range(1,self.n_procs):
                mm += mm_dict['m%d' % i_proc]
                bb += mm_dict['b%d' % i_proc]

        coeffs = np.linalg.solve(mm, bb)
        assert len(coeffs) == n_matrix

        fit_chisq = np.zeros(self.n_training, dtype=float)
        for i_proc in range(self.n_procs):
            i_start = self.pt_dexes[i_proc][0]
            i_end = self.pt_dexes[i_proc][1]
            pt_powers = self.pt_powers[i_proc]
            for i_matrix in range(n_matrix):
                fit_chisq[i_start:i_end] += coeffs[i_matrix]*pt_powers[i_matrix:len(pt_powers):self.dim*self.order]

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

        n_offenders = len(quad_4[0])+len(quad_2[0])
        n_non_offenders = len(quad_1[0])+len(quad_3[0])
        max_mismatch = max(chi_wrong[quad_4].max(), chi_wrong[quad_2].max())

        if max_mismatch<self.best_mismatch:
            self.best_mismatch = max_mismatch
            self.best_coeffs = coeffs

            with open('coeffs_test.txt', 'w') as out_file:
                out_file.write('# %d max_mismatch %e n_offenders %d\n' %
                               (self.i_iteration, max_mismatch, n_offenders))
                out_file.write('# dim %d order %d\n' % (self.dim, self.order))
                out_file.write('# i_matrix = i_dim*order + i_order\n')
                out_file.write('# chisq_min %e\n' % chisq_min)
                out_file.write('# center\n')
                out_file.write('# ')
                for cc in self.center:
                    out_file.write('%e ' % cc)
                out_file.write('\n')
                for cc in self.best_coeffs:
                    out_file.write('%e\n' % cc)
            with open('fit_chisq.txt', 'w') as out_file:
                out_file.write('# training_chisq-chisq_min fit_chisq rr training-fit sigma_sq\n')
                for i_pt in range(fitter.n_training):
                    out_file.write('%e %e %e %e %e\n' %
                                   (self.training_chisq[i_pt]-self.chisq_min,
                                    fit_chisq[i_pt],
                                    self.training_r[i_pt],
                                    self.training_chisq[i_pt]-self.chisq_min-fit_chisq[i_pt],
                                    sigma_sq[i_pt]))

        print('\n%d max_mis %.4e best %.4e -- %d %d -- %.2e %.2e %.2e' %
              (self.i_iteration, max_mismatch, self.best_mismatch,
               n_offenders, n_non_offenders,
               sigma_sq.min(),np.median(sigma_sq),sigma_sq.max()))

        print('')

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

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_proc', type=int, default=1)
    parser.add_argument('--n_train', type=float, default=0.66)
    parser.add_argument('--n_iter', type=int, default=100)
    parser.add_argument('--order', type=int, default=8)
    args = parser.parse_args()

    order = args.order
    dim = 33
    fitter = CartoonFitter(order,dim)
    basis_file = 'ellipse_bases.txt'
    fitter.read_bases(basis_file)
    data_file = 'planck_out_high_q2_culled_cartoon.txt'
    fitter.read_data(data_file, args.n_train, args.n_proc)

    sigma_sq = np.ones(len(fitter.training_chisq), dtype=float)
    sigma_sq = np.abs(fitter.training_chisq-fitter.chisq_min-fitter.delta_chisq)
    sigma_sq /= (0.1*fitter.delta_chisq)
    sigma_sq += 1.0
    sigma_sq = np.where(sigma_sq<10.0, sigma_sq, 10.0)

    #i_matrix = i_dim*order + i_order

    t_iter_start = time.time()
    for iteration in range(args.n_iter):
        t_last_iter = time.time()-t_iter_start
        print('\nlast iteration took %.2e min\n' % (t_last_iter/60.0))
        t_iter_start = time.time()

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
