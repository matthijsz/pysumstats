import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .qcplots import _assert_plot_defaults


def qqplot(pvector, fig=None, ax=None, filename=None, figsize=(8, 8), pointcolor='black', title=None, linecolor='red'):
    """Function to generate a QQ-plot.

    :param pvector: 1D-array of p-values
    :param fig: matplotlib.pyplot figure object to plot to (if not specified a new figure will be created)
    :param ax: matplotlib.pyplot axis to plot to (if not specified a new figure will be created)
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str.
    :param figsize: Figure size in inches (width, height)
    :type figsize: (float, float)
    :param pointcolor: Color to use for points
    :type pointcolor: str.
    :param title: Main figure title.
    :type title: str.
    :param linecolor: Color for line x=y
    :type linecolor: str.
    :return: None, or (fig, ax)

    """
    assert isinstance(pvector, pd.Series) or isinstance(pvector, np.ndarray), 'pvector should be pd.Series or np.ndarray'
    assert len(pvector.shape) == 1, 'pvector should be 1D-array'
    _assert_plot_defaults(fig=fig, ax=ax, filename=filename, pointcolor=pointcolor, linecolor=linecolor, title=title,
                          figsize=figsize)

    pvector = pvector[~np.isnan(pvector)]
    pvector = pvector[~((pvector >= 1) | (pvector <= 0))]
    pvector.sort()
    o = -(np.log10(pvector))
    e = -(np.log10((np.array(range(1, (len(o) + 1))) - .5) / len(o)))
    if (fig is None) and (ax is None):
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300*(int(figsize[0]/8)**2), facecolor='w')
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    maxi = max(np.amax(e), np.amax(o))
    ax.set_xlim([0, maxi + 0.1])
    ax.set_ylim([0, maxi + 0.1])
    ax.set_ylabel('Observed -log10(' + r'$p$' + ')')
    ax.set_xlabel('Expected -log10(' + r'$p$' + ')')
    ax.scatter(e, o, c=pointcolor)
    ax.plot((0, maxi), (0, maxi), linecolor)
    if title is not None:
        ax.set_title(title)
    if filename is not None:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        return fig, ax