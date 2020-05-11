import numpy as np
import matplotlib.pyplot as plt


def manhattan(dataframe, filename=None, sigp=5e-8, sigcolor='black', sugp=1e-5, sugcolor='black',
              pointcolor=['midnightblue', 'goldenrod'], size=1, highlight=[],
              highlightcolor=['orange'], title=None, rainbow=False):
    '''Create a Manhattan plot.

    :param dataframe: pd.Dataframe containing the following columns: ['rsid', 'chr', 'bp', 'p'], or :class:`pysumstats.SumStats`
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
    :param pointcolor: Color, or list of colors to cycle through for plotting SNPs
    :type pointcolor: list
    :param size: Relative figure size
    :type size: int or float
    :param highlight: list of SNPs to highlight
    :type highlight: list.
    :param highlightcolor: Color, or list of colors to cycle through for highlighting SNPs
    :type highlightcolor: list.
    :param title: Main figure title
    :type title: list.
    :param rainbow: Enble rainbow colors
    :type rainbow: bool.
    :return: None, or (fig, ax)

    '''
    df = dataframe[['rsid', 'chr', 'bp', 'p']].copy()
    df.columns = map(str.lower, df.columns)
    if rainbow:
        pointcolor = ['#FF0000', '#FF4000', '#FF8000', '#FFBF00', '#FFFF00', '#BFFF00', '#80FF00', '#40FF00', '#00FF00',
                      '#00FF40', '#00FF80', '#00FFBF', '#00FFFF', '#00BFFF', '#0080FF', '#0040FF', '#0000FF', '#4000FF',
                      '#8000FF', '#BF00FF', '#FF00FF', '#FF00BF']
        highlightcolor = list(reversed(pointcolor))
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
    fig, ax = plt.subplots(1, 1, figsize=(12 * float(size), 6 * float(size)), dpi=300*(float(size)**2), facecolor='w')
    colors = pointcolor
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='logp', color=colors[num % len(colors)],
                   edgecolor=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].mean()))
    if len(highlight) > 0:
        for num, (name, group) in enumerate(dfhighlight_grouped):
            group.plot(kind='scatter', x='ind', y='logp', color=highlightcolor[num % len(highlightcolor)],
                       edgecolor=highlightcolor[num % len(highlightcolor)], ax=ax, marker='^', s=45 * size)
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

