# Patch notes

##### 08-05-2020 (v0.1)

 - Adapted to be a package rather then a module.
 - Added `low_ram` argument to SumStats to read/write data to disk rather than RAM, in case of memory issues.  


# Description

A python module for working with GWAS summary statistics data in Python. <br/>
The main body of the code is designed to make it easy to read summary statistics, perform QC, merge summary statistics and perform meta-analysis.<br/>
Meta-analysis can be performed with .meta() with inverse-variance weighted or samplesize-weighted methods.<br/>
GWAMA as described in [Baselmans, et al. (2019)](https://www.nature.com/articles/s41588-018-0320-8) can be performed using the .gwama() function in merged summary statistics. <br/>
This module is fully compatible with my [Python-QQMan module](https://github.com/matthijsz/qqman).
Note merging with low_memory enabled is still highly experimental.

# Usage

See sumstats_usage.txt for example usage.

# Installation

This package was made for Python 3.7. Clone the package directly from this github, or install with 

`pip3 install pysumstats`
