import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from collections.abc import Iterable
from .qcplots import _assert_plot_defaults


def manhattan(dataframe, fig=None, ax=None, filename=None, sigp=5e-8, sigcolor='black', sugp=1e-5, sugcolor='black',
              pointcolors=['midnightblue', 'goldenrod'], figsize=(12, 6), highlight=[],
              highlightcolors=['orange'], title=None, rainbow=False):
    """Create a Manhattan plot.

    :param dataframe: pd.Dataframe containing the following columns: ['rsid', 'chr', 'bp', 'p'], or :class:`pysumstats.SumStats`
    :param fig: matplotlib.pyplot figure object to plot to (if not specified a new figure will be created)
    :param ax: matplotlib.pyplot axis to plot to (if not specified a new figure will be created)
    :param filename: Path to store the figure to (defaults to return fig, ax objects)
    :type filename: str.
    :param sigp: Where to plot significant line(set to a negative number to remove)
    :type sigp: float
    :param sigcolor: Color to use for significant line
    :type sigcolor: str
    :param sugp: Where to plot significant line (set to a negative number to remove)
    :type sugp: float
    :param sugcolor: Color to use for suggestive line
    :type sugcolor: str
    :param pointcolors: List of colors to cycle through for plotting SNPs
    :type pointcolors: list
    :param figsize: Figure size in inches (width, height)
    :type figsize: (float, float)
    :param highlight: list of SNPs to highlight
    :type highlight: list.
    :param highlightcolors: List of colors to cycle through for highlighting SNPs
    :type highlightcolors: list.
    :param title: Main figure title
    :type title: list.
    :param rainbow: Enble rainbow colors
    :type rainbow: bool.
    :return: None, or (fig, ax)
    """
    assert isinstance(dataframe, pd.DataFrame), 'dataframe should be a pandas.DataFrame object'
    assert np.all([x in dataframe.columns for x in ['chr', 'bp', 'p']]), "dataframe should contain ['chr', 'bp', 'p'] columns"
    assert isinstance(sigp, int) or isinstance(sigp, float), "sigp should be float"
    assert isinstance(sugp, int) or isinstance(sugp, float), "sigp should be float"
    assert is_color_like(sugcolor), "sugcolor should be interpretable as matplotlib color"
    assert is_color_like(sigcolor), "sigcolor should be interpretable as matplotlib color"
    assert isinstance(pointcolors, Iterable), "pointcolor should be list"
    for p in pointcolors:
        assert is_color_like(p), "values in highlightcolor should be interpretable as matplotlib color"
    assert isinstance(highlightcolors, Iterable), "highlightcolor should be list"
    for p in highlightcolors:
        assert is_color_like(p), "values in highlightcolor should be interpretable as matplotlib color"
    assert (title is None) or isinstance(title, str), "title should be None or str"
    assert isinstance(rainbow, bool), "rainbow should be bool"
    _assert_plot_defaults(fig=fig, ax=ax, filename=filename, title=title, figsize=figsize)

    if 'rsid' in dataframe.columns:
        df = dataframe[['rsid', 'chr', 'bp', 'p']].copy()
    else:
        df = dataframe[['chr', 'bp', 'p']].copy()
    df.columns = map(str.lower, df.columns)
    if rainbow:
        pointcolors = ['#FF0000', '#FF4000', '#FF8000', '#FFBF00', '#FFFF00', '#BFFF00', '#80FF00', '#40FF00', '#00FF00',
                      '#00FF40', '#00FF80', '#00FFBF', '#00FFFF', '#00BFFF', '#0080FF', '#0040FF', '#0000FF', '#4000FF',
                      '#8000FF', '#BF00FF', '#FF00FF', '#FF00BF']
        highlightcolors = list(reversed(pointcolors))
    df['logp'] = -(np.log10(df['p']))
    df['chromosome'] = df['chr'].astype('category')
    df['chromosome'] = df['chromosome'].cat.set_categories(['ch-%i' % i for i in range(23)], ordered=True)
    df.sort_values('chr', inplace=True)
    df['bpadd'] = 0
    add = 0
    for c in range(2, 24):
        add += df.loc[df['chr'] == (c - 1), 'bp'].max()
        df.loc[df['chr'] == c, 'bpadd'] = add
    df['ind'] = df['bp'] + df['bpadd']
    if len(highlight) > 0:
        dfhighlight = df[df['rsid'].isin(highlight)]
        dfhighlight_grouped = dfhighlight.groupby('chr')
    df_grouped = df.groupby('chr')
    if (fig is None) and (ax is None):
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=300*(float(figsize[0]/12)**2), facecolor='w')
    colors = pointcolors
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='logp', color=colors[num % len(colors)],
                   edgecolor=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].mean()))
    if len(highlight) > 0:
        for num, (name, group) in enumerate(dfhighlight_grouped):
            group.plot(kind='scatter', x='ind', y='logp', color=highlightcolors[num % len(highlightcolors)],
                       edgecolor=highlightcolors[num % len(highlightcolors)], ax=ax, marker='^', s=45)
    if sigp > 0:
        ax.plot((0, (df['ind'].max() * 1.005)), (-(np.log10(float(sigp))), -(np.log10(float(sigp)))), ls='--', lw=1.2, color=sigcolor)
    if sugp > 0:
        ax.plot((0, (df['ind'].max() * 1.005)), (-(np.log10(float(sugp))), -(np.log10(float(sugp)))), ls='--', lw=0.5, color=sugcolor)
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, (df['ind'].max() * 1.005)])
    ax.set_ylim([0, (np.amax(df['logp']) * 1.3)])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(' + r'$p$' + ')')
    if title is not None:
        ax.set_title(title)
    if filename is not None:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        return fig, ax

