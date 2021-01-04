"""
pysumstats, a python package for working with GWAS summary statistics

.. moduleauthor:: Matthijs D. van der Zee <m.d.vander.zee@vu.nl>
"""
__version__ = '0.5'
__author__ = 'Matthijs D. van der Zee'
__all__ = ['SumStats', 'cov_matrix_from_phenotype_file', 'plot', 'SumStatsWarning']

from . import plot
from .sumstats import (SumStats, MergedSumStats)
from .utils import cov_matrix_from_phenotype_file
from .exceptions import SumStatsWarning
