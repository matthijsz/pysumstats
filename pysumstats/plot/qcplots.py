import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import is_color_like
from collections.abc import Iterable
import matplotlib.lines as mlines
import scipy.stats as st
import warnings


def _assert_plot_defaults(twotailed=True, difference_cutoff=.1, differentcolor='red', differentlinecolor='red',
                          fig=None, ax=None, filename=None, pointcolor='black', linecolor='red', title=None,
                          figsize=(5, 5)):
    assert isinstance(twotailed, bool), "twotailed should be True or False"
    assert isinstance(difference_cutoff,
                      float) or difference_cutoff is None, "difference_cutoff should be float or None"
    assert is_color_like(differentcolor), "differentcolor should be interpretable as matplotlib color"
    assert is_color_like(differentlinecolor), "differentlinecolor should be interpretable as matplotlib color"
    assert (fig is None) or isinstance(fig, plt.Figure), "fig should be None or plt.Figure object"
    assert (ax is None) or isinstance(ax, plt.Axes), "ax should be None or plt.Axes object"
    assert (filename is None) or isinstance(filename, str), "filename should be None or str"
    assert is_color_like(pointcolor), "pointcolor should be interpretable as matplotlib color"
    assert is_color_like(linecolor), "linecolor should be interpretable as matplotlib color"
    assert (title is None) or isinstance(title, str), "title should be None or str"
    assert isinstance(figsize, Iterable), "figsize should be tuple"
    assert len(figsize) == 2, "figsize should be tuple of length 2"
    assert np.all([isinstance(x, float) for x in figsize]) or np.all(
        [isinstance(x, int) for x in figsize]), "figsize should contain floats"


def _seqdiffplot(ref, other, refname, othername, difference_cutoff, fig, ax, filename, pointcolor, differentcolor,
                 linecolor, differentlinecolor, title, figsize, xlim, ylim):
    """Base function for afplot, zzplot, pzplot

    :param ref: 1D-array of reference sequence
    :param other: 1D-array of other sequence
    :param refname: Name to use for reference sequence(x-axis label)
    :type refname: str
    :param othername: Name to use for other sequence (y-axis label)
    :type othername: str
    :param difference_cutoff: Cut-off to use for highlighting differences of sequence (to disable use None)
    :type difference_cutoff: None or float
    :param fig: matplotlib.pyplot figure object to plot to (if not specified a new figure will be created)
    :param ax: matplotlib.pyplot axis to plot to (if not specified a new figure will be created)
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str
    :param pointcolor: Color to use for points
    :type pointcolor: str
    :param differentcolor: Color to use for points that deviate given difference_cutoff
    :type differentcolor: str
    :param linecolor: Color to use for the line x=y
    :type linecolor: str
    :param differentlinecolor: Color to use for visualizing the difference_cutoff
    :type differentlinecolor: str
    :param title: Main figure title.
    :type title: str.
    :param figsize: Figure size
    :type figsize: (int, int)
    :return: None or (fig, ax)
    """
    if (fig is None) and (ax is None):
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    if len(ref) != len(other):
        raise ValueError('Lengths of input arrays do not match')
    if difference_cutoff is not None:
        dif = abs(ref - other)
        idx1 = np.where(dif < difference_cutoff)
        idx2 = np.where(dif >= difference_cutoff)
    else:
        idx1 = (np.array(list(range(ref.shape[0]))),)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_ylabel(othername)
    ax.set_xlabel(refname)
    ax.scatter(ref[idx1], other[idx1], c=pointcolor)
    if difference_cutoff is not None:
        ax.scatter(ref[idx2], other[idx2], c=differentcolor)
        l1 = mlines.Line2D([xlim[0], xlim[1] - difference_cutoff], [ylim[0] + difference_cutoff, ylim[1]],
                           color=differentlinecolor, linestyle='--')
        l2 = mlines.Line2D([xlim[0] + difference_cutoff, xlim[1]], [ylim[0], ylim[1] - difference_cutoff],
                           color=differentlinecolor, linestyle='--')
        ax.add_line(l1)
        ax.add_line(l2)
    l3 = mlines.Line2D(xlim, ylim, color=linecolor)
    ax.add_line(l3)
    if title is not None:
        ax.set_title(title)
    if filename is None:
        return fig, ax
    else:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        return None


def afplot(ref_af, other_af, refname='ref_EAF', othername='other_EAF', difference_cutoff=.1,
           fig=None, ax=None, filename=None, pointcolor='black', differentcolor='red', linecolor='black',
           differentlinecolor='red', title=None, figsize=(5, 5)):
    """Generate a plot of (differences in) allele frequencies of two summary statistics.

    :param ref_af: 1D-array of reference allele frequencies
    :param other_af: 1D-array of other allele frequencies
    :param refname: Name to use for reference allele frequencies (x-axis label)
    :type refname: str
    :param othername: Name to use for other allele frequencies (y-axis label)
    :type othername: str
    :param difference_cutoff: Cut-off to use for highlighting SNPs with different allele frequency (to disable use None)
    :type difference_cutoff: None or float
    :param fig: matplotlib.pyplot figure object to plot to (if not specified a new figure will be created)
    :param ax: matplotlib.pyplot axis to plot to (if not specified a new figure will be created)
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str
    :param pointcolor: Color to use for points
    :type pointcolor: str
    :param differentcolor: Color to use for points that deviate given difference_cutoff
    :type differentcolor: str
    :param linecolor: Color to use for the line x=y
    :type linecolor: str
    :param differentlinecolor: Color to use for visualizing the difference_cutoff
    :type differentlinecolor: str
    :param title: Main figure title.
    :type title: str.
    :param figsize: Figure size
    :type figsize: (int, int)
    :return: None or (fig, ax)
    """
    assert isinstance(ref_af, pd.Series) or isinstance(ref_af, np.ndarray), 'ref_af should be pd.Series or np.ndarray'
    assert len(ref_af.shape) == 1, 'ref_af should be 1D-array'
    assert isinstance(other_af, pd.Series) or isinstance(other_af, np.ndarray), 'other_af should be pd.Series or np.ndarray'
    assert len(other_af.shape) == 1, 'other_af should be 1D-array'
    assert len(ref_af) == len(other_af), 'lengths of ref_af and other_af should match'
    _assert_plot_defaults(difference_cutoff=difference_cutoff, differentcolor=differentcolor,
                          differentlinecolor=differentlinecolor, fig=fig, ax=ax, filename=filename,
                          pointcolor=pointcolor, linecolor=linecolor, title=title, figsize=figsize)

    xlim = [0, 1]
    ylim = [0, 1]
    return _seqdiffplot(ref_af, other_af, refname=refname, othername=othername, difference_cutoff=difference_cutoff,
                        fig=fig, ax=ax, filename=filename, pointcolor=pointcolor, differentcolor=differentcolor,
                        linecolor=linecolor, differentlinecolor=differentlinecolor, title=title, figsize=figsize, xlim=xlim,
                        ylim=ylim)


def pzplot(data, twotailed=True, difference_cutoff=.1, fig=None, ax=None, filename=None, pointcolor='black',
           differentcolor='red', linecolor='black', differentlinecolor='red', title=None, figsize=(5, 5)):
    """Generate a plot comparing the z-value as calculated from the p-value to the z-value as calculated from beta/se

    :param data: 2D-array containing the columns ['b', 'se', 'p']
    :param twotailed: True if p-value was computed from both ends of the distribution.
    :type twotailed: bool
    :param difference_cutoff: Cut-off to use for highlighting SNPs with different z-values (to disable use None)
    :type difference_cutoff: None or float
    :param fig: matplotlib.pyplot figure object to plot to (if not specified a new figure will be created)
    :param ax: matplotlib.pyplot axis to plot to (if not specified a new figure will be created)
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str
    :param pointcolor: Color to use for points
    :type pointcolor: str
    :param differentcolor: Color to use for points that deviate given difference_cutoff
    :type differentcolor: str
    :param linecolor: Color to use for the line x=y
    :type linecolor: str
    :param differentlinecolor: Color to use for visualizing the difference_cutoff
    :type differentlinecolor: str
    :param title: Main figure title.
    :type title: str.
    :param figsize: Figure size
    :type figsize: (int, int)
    :return: None or (fig, ax)
    """

    assert isinstance(data, pd.DataFrame), 'data should be pd.DataFrame'
    assert np.all([x in data.columns for x in ['b', 'se', 'p']]), "data should contain columns ['b', 'se', 'p']"
    _assert_plot_defaults(twotailed=twotailed, difference_cutoff=difference_cutoff, differentcolor=differentcolor,
                          differentlinecolor=differentlinecolor, fig=fig, ax=ax, filename=filename,
                          pointcolor=pointcolor, linecolor=linecolor, title=title, figsize=figsize)

    if twotailed:
        data['z_from_p'] = st.norm.ppf(1 - data['p'] / 2)
    else:
        data['z_from_p'] = st.norm.ppf(1 - data['p'])
    data['z_from_p'] *= np.sign(data['b'])
    data['z_from_b'] = data['b'] / data['se']
    data_ = data[['z_from_p', 'z_from_b']].copy()
    data.drop(axis=1, labels=['z_from_p', 'z_from_b'], inplace=True)
    data_['z_from_p'] = data_['z_from_p'].replace([np.inf, -np.inf], np.nan)
    data_['z_from_b'] = data_['z_from_b'].replace([np.inf, -np.inf], np.nan)
    data_ = data_.dropna()
    low = np.min([data_['z_from_b'].min(), data_['z_from_p'].min()])
    high = np.max([data_['z_from_b'].max(), data_['z_from_p'].max()])
    xlim = [low, high]
    ylim = [low, high]
    result = _seqdiffplot(data_['z_from_p'].values, data_['z_from_b'].values, refname=r'$Z$'+' from '+r'$p$',
                          othername=r'$Z = \frac{\beta}{se}$', difference_cutoff=difference_cutoff,
                          fig=fig, ax=ax, filename=filename, pointcolor=pointcolor, differentcolor=differentcolor,
                          linecolor=linecolor, differentlinecolor=differentlinecolor, title=title,
                          figsize=figsize, xlim=xlim, ylim=ylim)
    return result


def zzplot(data_x, data_y, xname='Z_x', yname='Z_y', twotailed=True, difference_cutoff=.5, fig=None, ax=None,
           filename=None, pointcolor='black', differentcolor='red', linecolor='black', differentlinecolor='red',
           title=None, figsize=(5, 5)):
    """Generate a plot comparing the z-values from two GWAS summary statistics

    :param data_x: 2D-array containing the column 'z', the columns ['b', 'se'] or the column 'p' (this is priority order)
    :param data_y: 2D-array containing the column 'z', the columns ['b', 'se'] or the column 'p' (this is priority order)
    :param twotailed: True if p-value was computed from both ends of the distribution.
    :type twotailed: bool
    :param xname: Name to use for z_values of data_x (x-axis label)
    :type xname: str
    :param yname: Name to use for z_values of data_y (y-axis label)
    :type yname: str
    :param difference_cutoff: Cut-off to use for highlighting SNPs with different z-values (to disable use None)
    :type difference_cutoff: None or float
    :param fig: matplotlib.pyplot figure object to plot to (if not specified a new figure will be created)
    :param ax: matplotlib.pyplot axis to plot to (if not specified a new figure will be created)
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str
    :param pointcolor: Color to use for points
    :type pointcolor: str
    :param differentcolor: Color to use for points that deviate given difference_cutoff
    :type differentcolor: str
    :param linecolor: Color to use for the line x=y
    :type linecolor: str
    :param differentlinecolor: Color to use for visualizing the difference_cutoff
    :type differentlinecolor: str
    :param title: Main figure title.
    :type title: str.
    :param figsize: Figure size
    :type figsize: (int, int)
    :return: None or (fig, ax)
    """

    assert isinstance(data_x, pd.DataFrame) or isinstance(data_x, pd.Series), 'data_x should be pd.DataFrame or pd.Series'
    if isinstance(data_x, pd.Series):
        assert data_x.name in ['z', 'p'], "data_x should contain columns 'z', ['b', 'se'] or 'p'"
    else:
        assert ('z' in data_x.columns) or ('p' in data_x.columns) or (np.all([x in data_x.columns for x in ['b', 'se']])), "data_x should contain columns 'z', ['b', 'se'] or 'p'"
    assert isinstance(data_y, pd.DataFrame) or isinstance(data_y,
                                                          pd.Series), 'data_x should be pd.DataFrame or pd.Series'
    if isinstance(data_y, pd.Series):
        assert data_y.name in ['z', 'p'], "data_x should contain columns 'z', ['b', 'se'] or 'p'"
    else:
        assert ('z' in data_y.columns) or ('p' in data_y.columns) or (
            np.all([x in data_y.columns for x in ['b', 'se']])), "data_x should contain columns 'z', ['b', 'se'] or 'p'"
    assert isinstance(xname, str), "xname should be str"
    assert isinstance(yname, str), "xname should be str"
    _assert_plot_defaults(twotailed=twotailed, difference_cutoff=difference_cutoff, differentcolor=differentcolor,
                          differentlinecolor=differentlinecolor, fig=fig, ax=ax, filename=filename,
                          pointcolor=pointcolor, linecolor=linecolor, title=title, figsize=figsize)


    data, zsource, keep = {}, {}, []
    for data_in, suffix in [[data_x, 'x'], [data_y, 'y']]:
        if 'z' in data_x.columns:
            data[suffix] = data_in['z'].values
            zsource[suffix] = 'z'
        elif ('b' in data_in.columns) and ('se' in data_in.columns):
            data[suffix] = (data_in['b'] / data_in['se']).values
            data[suffix] *= np.sign(data_in['b'])
            zsource[suffix] = 'b/se'
        elif 'p' in data_x.columns:
            if twotailed:
                data[suffix] = st.norm.ppf(1 - data_in['p'] / 2).values
            else:
                data[suffix] = st.norm.ppf(1 - data_in['p']).values
            zsource[suffix] = 'p'
        else:
            raise KeyError('No compatible columns found in data_{}'.format(suffix))
        keep = list(set(list(np.where(np.isfinite(data[suffix]))[0])))
    if zsource['x'] != zsource['y']:
        warnings.warn('Z value from data_x was obtained from {}, while Z value from data_y was obtained from {}'.format(
            zsource['x'], zsource['y']
        ))
    data['x'] = data['x'][keep]
    data['y'] = data['y'][keep]
    low = np.min([data['x'].min(), data['y'].min()])
    high = np.max([data['x'].max(), data['y'].max()])
    xlim = [low, high]
    ylim = [low, high]
    return _seqdiffplot(data['x'].values, data['y'].values, refname=xname,
                        othername=yname, difference_cutoff=difference_cutoff,
                        fig=fig, ax=ax, filename=filename, pointcolor=pointcolor, differentcolor=differentcolor,
                        linecolor=linecolor, differentlinecolor=differentlinecolor, title=title,
                        figsize=figsize, xlim=xlim, ylim=ylim)
