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
parser.add_argument('--sub_dir', type=str, default='T')
parser.add_argument('--x', type=int, default=None, nargs='+')
parser.add_argument('--d', type=int, default=None)
parser.add_argument('--c', type=float, default=21.03,
                    help='delta chisquared')
parser.add_argument('--limit', type=int, default=None)
parser.add_argument('--key_map', type=str, default=None)

import os

args = parser.parse_args()

if args.in_file is None:
    raise RuntimeError('must specify in_file')
if args.out_dir is None:
    raise RuntimeError('must specify out_dir')
if args.d is None:
    raise RutnimeError('must specify d')

use_sub_dir = True
if args.sub_dir.lower()[0] == 'f':
    use_sub_dir = False

if args.control == 'd12':
    args.control = os.path.join('/Users', 'danielsf', 'physics',
                                'MultiNest_v3.9', 'chains',
                                'gaussianJellyBean_d12_s99_n50000_t1.00e-03.txt')

if not os.path.isdir(args.out_dir):
    os.mkdir(args.out_dir)

if not isinstance(args.in_file, list):
    in_file_names = [args.in_file]
else:
    in_file_names = args.in_file

key_map_dict = {}
if args.key_map is not None:
    with open(args.key_map, 'r') as in_file:
        for line in in_file:
            params = line.strip().split()
            key_map_dict[int(params[0])] = params[1]

control_data = None
control_x = None
control_y = None

dalex_data = {}
for file_name in in_file_names:
    dalex_data[file_name] = None

for ix_dex in range(0,len(args.x),2):
    ix = args.x[ix_dex]
    iy = args.x[ix_dex+1]

    subdir = 'x%dy%d' % (ix, iy)
    if use_sub_dir:
        full_out_dir = os.path.join(args.out_dir, subdir)
        if not os.path.exists(full_out_dir):
            os.mkdir(full_out_dir)
    else:
        full_out_dir = args.out_dir

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
        alpha = 1.0
        if control_x is not None:
            plt.scatter(control_x, control_y, marker='o', color='k', s=5)
            alpha = 0.75
        plt.scatter(dx, dy, marker='x', color='r', s=5, alphat=alpha)

        if ix not in key_map_dict:
            plt.xlabel('$\\theta_{%d}$' % ix)
        else:
            plt.xlabel(key_map_dict[ix])

        if iy not in key_map_dict:
            plt.ylabel('$\\theta_{%d}$' % iy)
        else:
            plt.ylabel(key_map_dict[iy])

        plt.title('$\chi^2_{min} = %.2f$' % dmin)

        if args.limit is not None:
            npts=args.limit
        else:
            npts=len(dalex_data[file_name])

        out_name = os.path.join(full_out_dir,
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
    min_dex = np.argmin(good_data['chisq'])
    for ii in range(args.d):
        print("    dim %d min %e; range: %e -> %e" % (ii,
                                                      good_data['x%d' % ii][min_dex],
                                                      good_data['x%d' % ii].min(),
                                                      good_data['x%d' % ii].max()))
