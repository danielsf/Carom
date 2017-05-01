import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

nlive = 5000

import os

physics_dir = os.path.join('/Users', 'danielsf', 'physics')
multinest_file = os.path.join(physics_dir, 'MultiNest_v3.9', 'chains',
                              'nonGaussianLump_d12_s99_n5000_t1.00e-03.txt')

out_dir = os.path.join(physics_dir, 'Carom_drafts', 'figures')

from analyzeDalex import one_d_marginalized_posterior_from_multinest

data = None
i_fig = 0
for ii, ix in enumerate((0,3,6,9)):
    i_fig += 1
    if ii==0 or ii==6:
        plt.figsize = (30, 30)

    (xx, pp,
     data) = one_d_marginalized_posterior_from_multinest(multinest_file,
                                                         ix, 12,
                                                         data=data, dx=0.15)

    plt.subplot(2,2,i_fig)
    plt.plot(xx,pp)
    plt.xlabel('$\\theta_%d$' % ix)

plt.tight_layout()
plt.savefig(os.path.join(out_dir,'figure_9.eps'))
plt.close()
