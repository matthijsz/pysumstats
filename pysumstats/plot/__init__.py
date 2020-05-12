"""
pysumstats.plot, a subpackage for generating plots from summary statistics data

.. moduleauthor:: Matthijs D. van der Zee <m.d.vander.zee@vu.nl>
"""

__all__ = ['manhattan', 'qqplot', 'afplot', 'zzplot', 'pzplot']
__author__ = 'Matthijs D. van der Zee'

from .manhattan import manhattan
from .qq import qqplot
from .qcplots import (afplot, zzplot, pzplot)

