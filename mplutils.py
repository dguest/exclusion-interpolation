"""
Wrappers for basic matplotlib figures.
"""

import os

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.cm import get_cmap
from matplotlib.ticker import MultipleLocator

from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np

class Canvas:
    default_name = 'test.pdf'
    def __init__(self, out_path=None, figsize=(5.0,5.0*3/4), ext=None):
        self.fig = Figure(figsize)
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(1,1,1)
        self.out_path = out_path
        self.ext = ext

    def save(self, out_path=None, ext=None):
        output = out_path or self.out_path
        assert output, "an output file name is required"
        out_dir, out_file = os.path.split(output)
        if ext:
            out_file = '{}.{}'.format(out_file, ext.lstrip('.'))
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        self.canvas.print_figure(output, bbox_inches='tight')

    def __enter__(self):
        if not self.out_path:
            self.out_path = self.default_name
        return self
    def __exit__(self, extype, exval, extb):
        if extype:
            return None
        self.save(self.out_path, ext=self.ext)
        return True


# _________________________________________________________________________
# specific draw routines
_ax_size = 12
_text_size = 12

def draw2d_exclusion(can, z, axes, log=False, cb_label='', **kwargs):
    """
    Simple draw routine for 2d hist. Assumes an ndhist, and a
    canvas with attributes `ax` (the axes) and `fig` (the figure).
    """
    ax = can.ax
    fig = can.fig

    # colorbar tweaks
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="1.5%")

    xlims, ylims = axes[0].lims, axes[1].lims
    imextent = list(xlims) + list(ylims)
    cmap = get_cmap('hot')

    args = dict(aspect='auto', origin='lower', extent=imextent,
                cmap=cmap, interpolation='nearest')
    args.update(**kwargs)
    if log:
        args['norm'] = LogNorm()

    # set NaN values outside the normal range to prevent warnings
    vmin, vmax = np.nanmin(z), np.nanmax(z)
    z[np.isnan(z)] = vmin - 0.1*(vmax - vmin)

    im = ax.imshow(z, vmin=vmin, vmax=vmax, **args)
    cb = fig.colorbar(im, cax=cax)
    if cb_label:
        cb.set_label(cb_label, y=0.98, ha='right')
    set_axes(ax, axes)
    return im, cb

def set_axes(ax, axes, tick_mult=0.7):
    mlx = MultipleLocator(100)
    mly = MultipleLocator(20)
    xlims, ylims = axes[0].lims, axes[1].lims
    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)
    ax.get_xaxis().set_minor_locator(mlx)
    ax.get_yaxis().set_minor_locator(mly)
    ax.get_xaxis().set_ticks_position('both')
    ax.get_yaxis().set_ticks_position('both')
    ax.tick_params(labelsize=_ax_size, direction='in', which='both')
    ax.tick_params(which='minor', length=5*tick_mult)
    ax.tick_params(which='major', length=10*tick_mult)
    ax.set_xlabel(_ax_name(axes[0]), x=0.98, ha='right', size=_ax_size)
    ax.set_ylabel(_ax_name(axes[1]), y=0.98, ha='right', size=_ax_size)

# ________________________________________________________________________
# lower level draw utilities

def _ax_name(ax):
    nm, un = ax.name, ax.units
    return '{} [{}]'.format(nm, un) if un else nm

# bullshit utils

def helvetify():
    """
    Load 'Helvetica' default font (may be Arial for now)
    """
    from matplotlib import rc
    # 'Comic Sans MS', 'Trebuchet MS' works, Arial works in png...
    rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
