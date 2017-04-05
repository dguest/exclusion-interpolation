#!/usr/bin/env python3

"""
Draw exclusions from a text file
"""

from argparse import ArgumentParser
from mplutils import Canvas, set_axes, helvetify
import numpy as np
import os

def get_dict(textfile):
    out_dict = {}
    for line in textfile:
        mz, ma, xsec = line.split()
        out_dict[(int(mz), int(ma))] = float(xsec)
    return out_dict

def get_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('theory')
    parser.add_argument('expected')
    parser.add_argument('-z', '--mz-slice', type=int, default=1200, nargs='?')
    parser.add_argument('-e', '--ext', default='.pdf')
    parser.add_argument('-o', '--out-dir', default='plots')
    return parser.parse_args()

def get_axes():
    class Axis:
        def __init__(self, low, high, name, units):
            self.lims = (low, high)
            self.name = name
            self.units = units
        def get_pts(self):
            return np.arange(*self.lims, 1)

    return [Axis(300, 3000, r"$m_{Z'}$", "GeV"),
            Axis(300, 1000, r"$m_{A}$", "GeV")]

class ArbAx:
    def __init__(self, low, high):
        self.lims = (low, high)
        self.name = 'Theory / Expected'
        self.units = ''

# def get_xyz_arrays(grid_dict, key='observed'):
#     xyz = {(i[0], i[1], p[key]) for i, p in grid_dict.items()}
#     x, y, z = zip(*sorted(xyz))
#     return np.array(x), np.array(y), np.array(z)

def run():
    args = get_args()
    helvetify()
    with open(args.theory) as textfile:
        grid_theory = get_dict(textfile)
    with open(args.expected) as textfile:
        grid_xsec = get_dict(textfile)

    slice_points = {}
    for (mz, ma), val in grid_theory.items():
        if mz == args.mz_slice:
            slice_points[ma] = val / grid_xsec[(mz, ma)]

    masses, values = zip(*sorted(slice_points.items()))
    marr, varr = np.asarray(masses), np.asarray(values)

    mm = np.linspace(marr.min(), marr.max(), 1000)
    vv = np.interp(mm, marr, varr)
    log_vv = np.exp(np.interp(mm, marr, np.log(varr)))
    odir = args.out_dir
    if not os.path.isdir(odir):
        os.mkdir(odir)

    axes = get_axes()

    mzp = args.mz_slice
    with Canvas(f'{odir}/slice-{mzp}{args.ext}') as can:
        set_axes(can.ax, (axes[1], ArbAx(0, 2)))
        can.ax.plot(mm, vv, '-', color='blue', label='lin')
        can.ax.plot(mm, log_vv, '-', color='green', label='log')
        can.ax.axhline(1)
        can.ax.legend(framealpha=0, title=fr"$m_{{Z'}} = {mzp}$ GeV")
    with Canvas(f'{odir}/slice-{mzp}-log{args.ext}') as can:
        set_axes(can.ax, (axes[1], ArbAx(varr.min(), varr.max())))
        can.ax.plot(mm, vv, '-', color='blue', label='lin')
        can.ax.plot(mm, log_vv, '-', color='green', label='log')
        can.ax.axhline(1)
        can.ax.set_yscale('log')
        can.ax.legend(framealpha=0, title=fr"$m_{{Z'}} = {mzp}$ GeV")

if __name__ == '__main__':
    run()
