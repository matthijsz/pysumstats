'''
pysumstats.plot, a subpackage for generating plots from summary statistics data

.. moduleauthor:: Matthijs van der Zee <m.d.vander.zee@vu.nl>
'''
from .manhattan import manhattan
from .qq import qqplot

__all__ = ['manhattan', 'qqplot']