'''
pysumstats, a python package for working with GWAS summary statistics

.. moduleauthor:: Matthijs van der Zee <m.d.vander.zee@vu.nl>
'''
from .sumstats import (SumStats, MergedSumStats)
from .utils import cov_matrix_from_phenotype_file
from . import plot
__all__ = ['SumStats', 'cov_matrix_from_phenotype_file']
