import os
from analysisPythonScript import doAnalysis


if __name__ == "__main__":
    dim = 4
    delta_chisq = 9.5

    ct_list = range(10000, 40000+1, 5000)
    if len(ct_list)>9:
        raise RuntimeError("Cannot plot more than 9 ct steps: you have %d" % len(ct_list))

    ix_list = []
    iy_list = []
    for ii in range(dim):
        for jj in range(ii+1, dim):
            ix_list.append(ii)
            iy_list.append(jj)

    #ix_list = [0]
    #iy_list = [3]

    carom_dir = os.path.join('/Users', 'danielsf', 'physics', 'Carom')

    output_dir = os.path.join(carom_dir, 'figures', 'timeSeries151111')

    input_dir = os.path.join(carom_dir,'output')

    input_file = os.path.join(input_dir,'analysisTest','jellyBean_d4_s99_output.sav')

    control_names = {}
    for ix, iy in zip(ix_list, iy_list):
        control_file = os.path.join(carom_dir, 'controls', 'jellyBeanData')
        control_file = os.path.join(control_file, 'jellyBean_%d_%d_frequentistFullDrelative.txt' % (ix, iy))

        control_names['%d_%d' % (ix, iy)] = control_file

    doAnalysis(dim, delta_chisq, ix_list, iy_list, ct_list, input_file, control_names, output_dir)
