import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class _MultiWindowPlot:
    """Helper class for generating multi-windowed plots (mainly for MergedSumStats plots)

    :param n_plots: Number of windows to create
    :type n_plots: int
    :param nrows: Number of user-specified rows
    :type nrows: int
    :param ncols: Number of user-specified columns
    :type ncols: int
    :param figsize: User-specified figure size
    :type figsize: (int, int)
    :param shape: either 'rect' or 'square', used to aim for a rectangle, or square figure if user did not specify nrows and ncols
    :type shape: str
    """
    def __init__(self, n_plots, nrows=None, ncols=None, figsize=None, shape='rect'):
        if (nrows is None) and (ncols is None):
            if shape == 'rect':
                ncols = int(np.log2(n_plots/2))
                if ncols < 1:
                    ncols = 1
                nrows = int(np.ceil(n_plots/ncols))
            elif shape == 'square':
                nrows = int(np.sqrt(n_plots))
                if nrows < 1:
                    nrows = 1
                ncols = int(np.ceil(n_plots / nrows))
        elif nrows is None:
            nrows = int(np.ceil(n_plots / ncols))
        elif ncols is None:
            ncols = int(np.ceil(n_plots / nrows))
        else:
            if (nrows * ncols) < n_plots:
                raise ValueError('Not enough cells to plot to.')
        if figsize is None:
            if shape == 'rect':
                figsize = (ncols*8, nrows*4)
            elif shape == 'square':
                figsize = (ncols*5, nrows*5)
        self.nrows, self.ncols = nrows, ncols
        self.fig, self.axs = plt.subplots(nrows, ncols, figsize=figsize)
        self.fig.tight_layout(pad=2, h_pad=2.2)
        self.n_plot = 0

    def get_next_ax(self, byrow=True):
        """Get the next axis to plot to.

        :param byrow: Fill subplots by rows, if False fill subplots by column
        :type byrow: bool
        :return: (fig, ax)
        """
        if (self.nrows > 1) and (self.ncols > 1):
            if byrow:
                ax = self.axs[(self.n_plot % self.nrows), (self.n_plot // self.nrows)]
            else:
                ax = self.axs[(self.n_plot // self.ncols), (self.n_plot % self.ncols)]
        elif (self.nrows > 1) or (self.ncols > 1):
            ax = self.axs[self.n_plot]
        else:
            ax = self.axs
        self.n_plot += 1
        return self.fig, ax

    def finish(self, filename, dpi):
        """Finish the plot with given filename and dpi

        :param filename: None, or path to target file (if None, (fig, ax) will be returned
        :param dpi: dpi setting to use when path is specified
        :return: (fig, ax) or None
        """
        if filename is None:
            return self.fig, self.axs
        else:
            plt.savefig(filename, dpi=dpi)
            plt.close()
            return None


def cov_matrix_from_phenotype_file(dataframe, phenotypes=None):
    """Function to generate a covariance (cov_Z) matrix from a phenotype file.

    :param dataframe: pd.Dataframe containing phenotypic data
    :type dataframe: pd.Dataframe
    :param phenotypes: list of phenotypes to include
    :type phenotypes: list.
    :return: pd.Dataframe of covariance matrix.

    """
    if phenotypes is None:
        phenotypes = list(dataframe.columns)
    cov_matrix = pd.DataFrame(0, index=phenotypes, columns=phenotypes)
    for p1 in phenotypes:
        cov_matrix.loc[p1, p1] = 1
        n1 = len(dataframe.loc[dataframe[p1].notna(), :])
        for p2 in [x for x in phenotypes if x != p1]:
            ns = len(dataframe.loc[(dataframe[p1].notna() & dataframe[p2].notna()), :])
            n2 = len(dataframe.loc[dataframe[p2].notna(), :])
            r = dataframe[[p1, p2]].corr().iloc[0, 1]
            cov_matrix.loc[p1, p2] = (r * ns) / np.sqrt((n1 * n2))
    return cov_matrix
