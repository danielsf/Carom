import os
from analysisPythonScript import doAnalysis


if __name__ == "__main__":
    dim = 6
    delta_chisq = 12.6

    ct_list = range(10000, 40000+1, 5000)
    if len(ct_list)>9:
        raise RuntimeError("Cannot plot more than 9 ct steps: you have %d" % len(ct_list))

    ix_list = []
    iy_list = []
    for ii in range(dim):
        for jj in range(ii+1, dim):
            if ii!=3 and jj!=3:
                ix_list.append(ii)
                iy_list.append(jj)

    #ix_list = [0]
    #iy_list = [3]

    carom_dir = os.path.join('/Users', 'danielsf', 'physics', 'Carom')

    output_dir = os.path.join(carom_dir, 'figures', 'wmapScratch')

    input_dir = os.path.join(carom_dir,'output')

    input_file = os.path.join(input_dir,'scratch', 'wmap7_s122_output.sav')

    control_names = {}
    for ix, iy in zip(ix_list, iy_list):
        control_file = os.path.join(carom_dir, 'processedChains', 'wmap7')
        control_file = os.path.join(control_file, 'mcmc_test_150530_carom_30k_%d_%d_frequentist.sav' % (ix, iy))

        control_names['%d_%d' % (ix, iy)] = control_file

    doAnalysis(dim, delta_chisq, ix_list, iy_list, ct_list, input_file, control_names, output_dir, scatter_control=True)
