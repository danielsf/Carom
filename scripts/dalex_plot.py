import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse

from analyzeDalex import scatter_from_multinest_projection
from analyzeDalex import scatter_from_dalex

parser = argparse.ArgumentParser()
parser.add_argument('--in_file', type=str, default=None,
                    nargs='+')
parser.add_argument('--in_dir', type=str, default=None)
parser.add_argument('--control', type=str, default=None)
parser.add_argument('--out_dir', type=str, default=None)
parser.add_argument('--x', type=int, default=None, nargs='+')
parser.add_argument('--d', type=int, default=None)
parser.add_argument('--c', type=float, default=21.03,
                    help='delta chisquared')
parser.add_argument('--limit', type=int, default=None)

import os

args = parser.parse_args()

if args.in_file is None:
    raise RuntimeError('must specify in_file')
if args.out_dir is None:
    raise RuntimeError('must specify out_dir')
if args.d is None:
    raise RutnimeError('must specify d')

if not os.path.isdir(args.out_dir):
    os.mkdir(args.out_dir)

if not isinstance(args.in_file, list):
    in_file_names = [args.in_file]
else:
    in_file_names = args.in_file

control_data = None
control_x = None
control_y = None

dalex_data = {}
for file_name in in_file_names:
    dalex_data[file_name] = None

for ix_dex in range(0,len(args.x),2):
    ix = args.x[ix_dex]
    iy = args.x[ix_dex+1]

    if args.control is not None:
        (control_x,
         control_y,
         control_data) = scatter_from_multinest_projection(args.control,
                                                           args.d,
                                                           ix, iy,
                                                           data=control_data)

    for file_name in in_file_names:
        if args.in_dir is not None:
            full_name = os.path.join(args.in_dir, file_name)
        else:
            full_name = file_name

        (dx, dy,
         dmin, dtarget,
         ddata) = scatter_from_dalex(full_name, args.d, ix, iy,
                                     delta_chi=args.c,
                                     data=dalex_data[file_name],
                                     limit=args.limit)

        dalex_data[file_name] = ddata

        plt.figsize = (30,30)
        if control_x is not None:
            plt.scatter(control_x, control_y, marker='o', color='k', s=10)
        plt.scatter(dx, dy, marker='x', color='r', s=10)
        plt.xlabel('$\\theta_{%d}$' % ix)
        plt.ylabel('$\\theta_{%d}$' % iy)

        if args.limit is not None:
            npts=args.limit
        else:
            npts=len(dalex_data[file_name])

        out_name = os.path.join(args.out_dir,
                                '%s_%d_%d_%.2e.png' % (file_name, ix, iy, npts))

        plt.savefig(out_name)
        plt.close()

import numpy as np
print('npts %d' % npts)
for file_name in dalex_data:
    data_arr = dalex_data[file_name][:npts]
    good_dexes = np.where(data_arr['chisq']<=dtarget)
    good_data = data_arr[good_dexes]
    print("stats for %s" % file_name)
    print("    chi^2 min: %e" % good_data['chisq'].min())
    for ii in range(args.d):
        print("    range %d: %e -> %e" % (ii,
                                          good_data['x%d' % ii].min(),
                                          good_data['x%d' % ii].max()))
