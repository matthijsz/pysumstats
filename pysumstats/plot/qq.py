import numpy as np
import matplotlib.pyplot as plt


def qqplot(pvector, filename=None, size=1, pointcolor='black', title=None, linecolor='red'):
    """Function to generate a QQ-plot.

    :param pvector: 1D-array of p-values
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str.
    :param size: Relative figure size (default=1)
    :type size: int or float
    :param pointcolor: Color to use for points
    :type pointcolor: str.
    :param title: Main figure title.
    :type title: str.
    :param linecolor: Color for line x=y
    :type linecolor: str.
    :return: None, or (fig, ax)

    """
    pvector = pvector[~np.isnan(pvector)]
    pvector = pvector[~((pvector >= 1) | (pvector <= 0))]
    pvector.sort()
    o = -(np.log10(pvector))
    e = -(np.log10((np.array(range(1, (len(o) + 1))) - .5) / len(o)))
    fig, ax = plt.subplots(1, 1, figsize=(8 * int(size), 8 * int(size)), dpi=300*(int(size)**2), facecolor='w')
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