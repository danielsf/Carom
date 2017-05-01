from __future__ import with_statement
import argparse
import os

from analyzeDalex import scatter_from_multinest_projection

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, default=None)
    parser.add_argument('--outfile', type=str, default=None)
    parser.add_argument('--clobber', type=bool, default=False)
    parser.add_argument('--dim', type=int, default=None)
    parser.add_argument('--x', type=int, nargs='+', default=None)

    args = parser.parse_args()
    if args.infile is None or args.outfile is None:
        raise RuntimeError("input %s output %s" % (args.infile, args.outfile))

    if args.dim is None or args.x is None:
        raise RuntimeError('dim %s x %s' % (args.dim, args.x))

    if not args.clobber and os.path.exists(args.outfile):
        raise RuntimeError("Cannot write %s; already exists" % args.outfile)

    mx, my, data = scatter_from_multinest_projection(args.infile,
                                                     args.dim,
                                                     args.x[0], args.x[1])


    with open(args.outfile, 'w') as output_file:
        for xx, yy in zip(mx, my):
            output_file.write('%e %e\n' % (xx, yy))

    for ix in range(args.dim):
        print '%d %e %e' % (ix, data['x%d'%ix].min(), data['x%d'%ix].max())
