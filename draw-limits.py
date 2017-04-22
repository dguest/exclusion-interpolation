#!/usr/bin/env python3

"""
Draw exclusions from a text file
"""

from argparse import ArgumentParser
from mplutils import Canvas, draw2d_exclusion, set_axes, helvetify
import numpy as np
import os

from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.stats import norm
from scipy.spatial import Delaunay

def get_grid_dict(textfile):
    file_iter = iter(textfile)
    titles = next(file_iter).split()
    sig_keys = titles[1:]
    grid_dict = {}
    for line in file_iter:
        line_vals = {n: v for n, v in zip(titles, line.split())}
        model, mzs, mas = line_vals['file'].split('_')
        mz, ma = int(mzs), int(mas)
        grid_dict[(mz, ma)] = {n: float(line_vals[n]) for n in sig_keys}
    return grid_dict

def get_xsec_dict(textfile):
    out_dict = {}
    for line in textfile:
        mz, ma, xsec = line.split()
        out_dict[(int(mz), int(ma))] = {'xsec': float(xsec)}
    return out_dict

def get_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('theory')
    parser.add_argument('expected')
    parser.add_argument('ex_u1')
    parser.add_argument('ex_d1')
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

def get_xyz_arrays(grid_dict, key='observed'):
    xyz = {(i[0], i[1], p[key]) for i, p in grid_dict.items()}
    x, y, z = zip(*sorted(xyz))
    return np.array(x), np.array(y), np.array(z)

def run():
    args = get_args()
    helvetify()
    with open(args.theory) as textfile:
        grid_theory = get_xsec_dict(textfile)
    x, y, z = get_xyz_arrays(grid_theory, 'xsec')
    bdict = {}
    bands = [('exp', args.expected),
             ('up', args.ex_u1),
             ('dn', args.ex_d1)]
    for name, band in bands:
        with open(band) as textfile:
            grid_xsec = get_xsec_dict(textfile)
        bdict[name] = get_xyz_arrays(grid_xsec, 'xsec')[2]

    xax, yax = get_axes()
    axes = [xax, yax]
    x_pts, y_pts = xax.get_pts(), yax.get_pts()
    xx, yy = np.meshgrid(x_pts, y_pts)
    grid = get_grid(x, y)
    th_dict = {
        'lin': _interpolate_linear(grid, z, xx, yy),
        'log': _interpolate_log(grid, z, xx, yy)
    }
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    odir = args.out_dir
    with Canvas(f'{odir}/triangles{args.ext}') as can:
        set_axes(can.ax, axes)
        can.ax.triplot(x, y, grid.simplices.copy())

    cbx = 'Cross Section [pb]'
    cb_th = cbx + ' (theory)'
    with Canvas(f'{odir}/theory{args.ext}') as can:
        zz = th_dict['lin']
        draw2d_exclusion(can, zz, [xax, yax], log=True, cb_label=cb_th)
                         # vmin=zz.min(), vmax=zz.max(), log=True)
        can.ax.plot(x, y, '.')
    with Canvas(f'{odir}/theory-log{args.ext}') as can:
        zz = th_dict['log']
        draw2d_exclusion(can, zz, [xax, yax], log=True, cb_label=cb_th)
                         # vmin=zz.min(), vmax=zz.max(), log=True)
        can.ax.plot(x, y, '.')

    z_grids = {}
    for name, z_vals in bdict.items():
        lin = _interpolate_linear(grid, z_vals, xx, yy)
        log = _interpolate_log(grid, z_vals, xx, yy)
        z_grids[name] = {'lin': lin, 'log': log}

        if name != 'exp': continue
        cb_exp = cbx + ' (expected)'
        with Canvas(f'{odir}/{name}{args.ext}') as can:
            draw2d_exclusion(can, lin, [xax, yax], log=True, cb_label=cb_exp)
            can.ax.plot(x, y, '.')
        with Canvas(f'{odir}/{name}-log{args.ext}') as can:
            draw2d_exclusion(can, log, [xax, yax], log=True, cb_label=cb_exp)
            can.ax.plot(x, y, '.')

        ratio_lab = (r'$\frac{\mathrm{Theory} - \mathrm{Expected}}'
                     r'{\mathrm{Theory} + \mathrm{Expected}}$')
        with Canvas(f'{odir}/theory-minus-{name}{args.ext}') as can:
            lab = ratio_lab + ' (lin interp)'
            zz = (th_dict['lin'] - lin) / (th_dict['lin'] + lin)
            im, cb = draw2d_exclusion(can, zz, [xax, yax],
                                      log=False, cb_label=lab)
            cs = can.ax.contour(xx, yy, th_dict['lin'] - lin, [0])
            cb.add_lines(cs)
            can.ax.plot(x, y, '.')
        with Canvas(f'{odir}/theory-minus-{name}-log{args.ext}') as can:
            lab = ratio_lab + ' (log interp)'
            zz = (th_dict['log'] - log) / (th_dict['log'] + log)
            im, cb = draw2d_exclusion(can, zz, [xax, yax],
                                      log=False, cb_label=lab)
            cs = can.ax.contour(xx, yy, th_dict['log'] - log, [0])
            cb.add_lines(cs)
            can.ax.plot(x, y, '.')

    for scale in ['lin','log']:
        zz_th = th_dict[scale]

        lowp = z_grids['dn'][scale] - zz_th
        highp = z_grids['up'][scale] - zz_th
        nom = z_grids['exp'][scale] - zz_th
        with Canvas(f'{odir}/bands-{scale}{args.ext}') as can:
            set_axes(can.ax, axes)
            zp = np.maximum( (lowp - 0), -(highp - 0))
            can.ax.contourf(xx, yy, zp, [-1, 0],
                            colors=['lime'], zorder=0)
            can.ax.contour(xx, yy, nom, [0], colors=['k'],
                           linestyles=['--'])
            can.ax.plot(x, y, '.')

        with Canvas(f'{odir}/exp-minus-theory-{scale}{args.ext}') as can:
            set_axes(can.ax, axes)
            for name, var in z_grids.items():
                style = '-' if name == 'exp' else ':'
                zz = var[scale] - zz_th
                can.ax.contour(xx, yy, zz, [0], colors=['k'],
                               linestyles=[style])
                can.ax.plot(x, y, '.')


def _interpolate_linear(pts, z, xp, yp):
    lin = LinearNDInterpolator(pts, z)
    interp_points = np.vstack((xp.flatten(), yp.flatten())).T
    zp = lin(interp_points).reshape(xp.shape)
    return zp

def _interpolate_normal(xy_grid, z, xp, yp):
    lin = LinearNDInterpolator(xy_grid, norm.ppf(z))
    interp_points = np.vstack((xp.flatten(), yp.flatten())).T
    interp_z = lin(interp_points).reshape(xp.shape)
    # avoid passing nan into the norm.cdf function
    nans = np.isnan(interp_z)
    zp = np.empty(interp_z.shape)
    zp[nans] = np.nan
    zp[~nans] = norm.cdf(interp_z[~nans])
    return zp

def get_grid(x, y):
    pts = np.vstack((x - 1e-9*y,y)).T
    return Delaunay(pts)

def _interpolate_log(xy_grid, z, xp, yp):
    lin = LinearNDInterpolator(xy_grid, np.log(z))
    interp_points = np.vstack((xp.flatten(), yp.flatten())).T
    interp_z = lin(interp_points).reshape(xp.shape)
    return np.exp(interp_z)

if __name__ == '__main__':
    run()
