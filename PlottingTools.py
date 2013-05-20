from __future__ import division
__author__ = 'will'
import pandas
import numpy as np
from matplotlib import pyplot as plt
from pylab import get_cmap
from types import StringType


def make_heatmap(data, col_labels, row_labels, colormap=None, **kwargs):

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
    plt.hold(True)
    return fig


def make_heatmap_df(dataframe, **kwargs):

    row_labels = dataframe.index
    col_labels = dataframe.columns
    return make_heatmap(dataframe.values, col_labels, row_labels, **kwargs)
