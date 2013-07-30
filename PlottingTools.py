from __future__ import division
__author__ = 'will'
import pandas
import numpy as np
from matplotlib import pyplot as plt
from pylab import get_cmap
from types import StringType


def make_heatmap(data, col_labels, row_labels, make_grid=True,
                 grid_kwargs=None, colormap=None, **kwargs):

    if make_grid and (grid_kwargs is None):
        grid_kwargs = {'color': 'w', 'lw': 2}

    fig = plt.figure(**kwargs)
    plt.hold(True)
    if colormap is None:
        cmap = None
    elif type(colormap) == StringType:
        cmap = get_cmap(colormap)
    else:
        cmap = colormap

    plt.imshow(data, interpolation='nearest', cmap=cmap, aspect='auto')

    plt.yticks(range(len(row_labels)), row_labels)
    plt.xticks(range(len(col_labels)), col_labels, rotation=90)

    if grid_kwargs:
        xpos = np.arange(len(col_labels))+0.5
        ypos = np.arange(len(row_labels))+0.5
        add_grid(plt.gca(), xpos, ypos, **grid_kwargs)

    plt.hold(True)
    return fig


def pos_to_rectangle(xmin, xmax, ymin, ymax):
    """Converts (Xmin, Xmax) (Ymin, Ymax) pairs into
    [left, bottom, width, height] needed by matplotlib rectangles.
    """

    left = xmin
    bottom = ymin
    width = xmax - xmin
    height = ymax - ymin
    return [left, bottom, width, height]


def generate_broken_axis(major_xlim, minor_xlim, major_ylim, minor_ylim,
                         ratio=0.3, fig_specs={}, axes_specs={},
                         edge_spacing=0.05, broken_spacing=0.05):

    if all(x == y for x, y in zip(major_xlim, minor_xlim)):
        is_vertical = True
        minor_above = minor_ylim[0] > major_ylim[0]
    elif all(x == y for x, y in zip(major_ylim, minor_ylim)):
        is_vertical = False
        minor_above = minor_xlim[0] > major_xlim[0]
    else:
        raise ValueError('limits must be the same on at least one axis!')


    fig = plt.figure(**fig_specs)

    es = edge_spacing
    bs = broken_spacing
    #assumes a 0.05 spacing between axes and edges
    # edge_spacing+x+broken_spacing+y+edge_spacing=1
    # where x and y are axes-sizes
    # y = ratio*x since the smaller axis is a ratio of the larger
    major_size = (1 - 2*edge_spacing - broken_spacing)/(ratio+1)
    minor_size = major_size*ratio

    #shape = [left, bottom, width, height]
    if is_vertical and minor_above:
        major_shape = pos_to_rectangle(es, 1-es, es, major_size + es)
        minor_shape = pos_to_rectangle(es, 1-es, major_size + es + bs, 1 - es)
    elif is_vertical and not minor_above:
        minor_shape = pos_to_rectangle(es, 1-es, 0.05, 0.05 + minor_size)
        major_shape = pos_to_rectangle(es, 1-es, minor_size + 0.1, 0.95)
    elif not is_vertical and minor_above:
        major_shape = pos_to_rectangle(0.05, 0.05 + major_size, 0.05, 0.95,)
        minor_shape = pos_to_rectangle(0.1 + major_size, 0.95, 0.05, 0.95,)
    elif not is_vertical and not minor_above:
        minor_shape = pos_to_rectangle(0.05, 0.05 + minor_size, 0.05, 0.95)
        major_shape = pos_to_rectangle(minor_size + 0.1, 0.95, 0.05, 0.95)


def add_grid(ax, xgrid_pos, ygrid_pos, **kwargs):
    """Adds a 'grid' to the provided axis. Any kwargs are passed as line-specs."""

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    ax.hlines(ygrid_pos, xmin, xmax, **kwargs)
    ax.vlines(xgrid_pos, ymin, ymax, **kwargs)


def make_heatmap_df(dataframe, **kwargs):

    row_labels = dataframe.index
    col_labels = dataframe.columns
    return make_heatmap(dataframe.values, col_labels, row_labels, **kwargs)
