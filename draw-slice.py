#!/usr/bin/env python3

"""
Draw exclusions from a text file
"""

from argparse import ArgumentParser
from mplutils import Canvas, set_axes, helvetify
import numpy as np
from scipy.interpolate import interp1d
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
    parser.add_argument('-s', '--do-splines', action='store_true')
    parser.add_argument('--hbb-br', nargs='?', default=0.571)
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
class XSecAx:
    def __init__(self, low, high):
        self.lims = (low, high)
        self.name = 'Cross Section'
        self.units = 'pb'

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

    th_points = {}
    ex_points = {}
    for (mz, ma), val in grid_theory.items():
        if mz == args.mz_slice:
            th_points[ma] = val
            ex_points[ma] = grid_xsec[(mz, ma)]

    masses, th_values = zip(*sorted(th_points.items()))
    masses, ex_values = zip(*sorted(ex_points.items()))
    marr = np.array(masses)
    tarr = np.array(th_values) * args.hbb_br
    earr = np.array(ex_values)
    varr = tarr / earr

    mm = np.linspace(marr.min(), marr.max(), 1000)
    def logint(vals):
        return np.exp(np.interp(mm, marr, np.log(vals)))
    def spl_int(vals):
        return interp1d(marr, vals, kind='quadratic')(mm)
    vv = np.interp(mm, marr, varr)
    log_vv = logint(varr)
    tt = np.interp(mm, marr, tarr)
    log_tt = logint(tarr)
    spl_tt = spl_int(tarr)
    ee = np.interp(mm, marr, earr)
    log_ee = logint(earr)
    spl_ee = spl_int(earr)
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

    with Canvas(f'{odir}/slice-xsec-{mzp}{args.ext}') as can:
        set_axes(can.ax, (axes[1], XSecAx(tt.min(), tt.max())))
        can.ax.plot(mm, tt, '-', color='red', label='Theory')
        can.ax.plot(mm, ee, '-', color='green', label='Expected')
        can.ax.plot(mm, log_tt, '--', color='red', label='Theory Log')
        can.ax.plot(mm, log_ee, '--', color='green', label='Expected Log')
        if args.do_splines:
            can.ax.plot(mm, spl_tt, ':', color='red', label='Theory Spl')
            can.ax.plot(mm, spl_ee, ':', color='green', label='Expected Spl')
        can.ax.axhline(1)
        can.ax.legend(framealpha=0, title=fr"$m_{{Z'}} = {mzp}$ GeV")
    with Canvas(f'{odir}/slice-xsec-{mzp}-log{args.ext}') as can:
        set_axes(can.ax, (axes[1], XSecAx(tt.min(), tt.max())))
        can.ax.plot(mm, tt, '-', color='red', label='Theory')
        can.ax.plot(mm, ee, '-', color='green', label='Expected')
        can.ax.plot(mm, log_tt, '--', color='red', label='Theory Log')
        can.ax.plot(mm, log_ee, '--', color='green', label='Expected Log')
        if args.do_splines:
            can.ax.plot(mm, spl_tt, ':', color='red', label='Theory Spl')
            can.ax.plot(mm, spl_ee, ':', color='green', label='Expected Spl')
        can.ax.axhline(1)
        can.ax.set_yscale('log')
        can.ax.legend(framealpha=0, title=fr"$m_{{Z'}} = {mzp}$ GeV")

if __name__ == '__main__':
    run()
