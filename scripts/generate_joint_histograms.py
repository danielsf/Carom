import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from analyzeDalex import make_histogram

physics_dir = os.path.join("/Users", "danielsf", "physics")
dalex_dir = os.path.join(physics_dir, "Carom", "output",
                            "workspace")

multinest_dir = os.path.join(physics_dir, "MultiNest_v3.9", "chains")
fig_dir = os.path.join(physics_dir, "Carom_drafts", "figures")

multinest_nlive = 300
dalex_seed = 90

dchi = 0.5
chi_cutoff = 250.0
chi_min = 94.0

multinest_name = os.path.join(multinest_dir,
                'gaussianJellyBean_d12_s99_n%d_t1.00e-03_carom.sav'
                % multinest_nlive)

dalex_name = os.path.join(dalex_dir, 'jellyBean_d12_s%d_output.sav' % dalex_seed)

dtype_list = []
for ii in range(12):
    dtype_list.append(("p%d" % ii, float))
dtype_list.append(("chisq", float))
dtype_list.append(("j1", int))
dtype_list.append(("j2", int))
dtype_list.append(("log", int))
dtype = np.dtype(dtype_list)

dalex_data = np.genfromtxt(dalex_name, dtype=dtype)
multinest_data = np.genfromtxt(multinest_name, dtype=dtype)

header_list = []
label_list = []

plt.figsize = (30, 30)
limit = 200000
dx_200k, dy_200k = make_histogram(dalex_data['chisq'][:limit],
                                  dchi, chi_cutoff,
                                  cumulative=False,
                                  min_val=chi_min)

dy_200k = dy_200k/dchi

hh, = plt.plot(dx_200k, dy_200k)
header_list.append(hh)
label_list.append('Dale$\chi$: %.2e $\chi^2$ calls; $\chi^2_{min}=%.2f$'
                  % (limit, dalex_data['chisq'][:limit].min()))


dx_max, dy_max = make_histogram(dalex_data['chisq'],
                                dchi, chi_cutoff,
                                cumulative=False,
                                min_val=chi_min)

dy_max = dy_max/dchi

hh, = plt.plot(dx_max, dy_max)
header_list.append(hh)
label_list.append('Dale$\chi$: %.2e $\chi^2$ calls; $\chi^2_{min}=%.2f$'
                  % (len(dalex_data['chisq']), dalex_data['chisq'].min()))

mx_200k, my_200k = make_histogram(multinest_data['chisq'][:limit],
                                  dchi, chi_cutoff,
                                  cumulative=False,
                                  min_val=chi_min)

my_200k = my_200k/dchi

hh, = plt.plot(mx_200k, my_200k)
header_list.append(hh)
label_list.append('MultiNest: %.2e $\chi^2$ calls; $\chi^2_{min}=%.2f$'
                  % (limit, multinest_data['chisq'][:limit].min()))


mx_max, my_max = make_histogram(multinest_data['chisq'],
                                dchi, chi_cutoff,
                                cumulative=False,
                                min_val=chi_min)

my_max = my_max/dchi

hh, = plt.plot(mx_max, my_max)
header_list.append(hh)
label_list.append('MultiNest: %.2e $\chi^2$ calls; $\chi^2_{min}=%.2f$'
                  % (len(multinest_data['chisq']),
                     multinest_data['chisq'].min()))

plt.legend(header_list, label_list, fontsize=10)
plt.xlabel('$\chi^2$')
plt.ylabel('dN/d$\chi^2$')
plt.xlim(dalex_data['chisq'].min(), chi_cutoff)

plt.savefig(os.path.join(fig_dir, 'figure_13.png'))
