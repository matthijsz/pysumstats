[![Documentation Status](https://readthedocs.org/projects/pysumstats/badge/?version=latest)](https://pysumstats.readthedocs.io/en/latest/?badge=latest)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![PyPI version](https://badge.fury.io/py/pysumstats.svg)](https://badge.fury.io/py/pysumstats)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.org/matthijsz/pysumstats.svg?branch=master)](https://travis-ci.org/matthijsz/pysumstats)

# Patch notes

##### 30-11-2020 (v0.4.5)
 - Fix: Merge using a list of SumStats now propagates arguments properly
 - Fix: Meta analyzing now properly sums N across all summary statistics
 

#### Previous
Older patchnodes can be found in [PATCHNOTES.md](PATCHNOTES.md)

# Description

A python package for working with GWAS summary statistics data in Python. <br/>
This package is designed to make it easy to read summary statistics, perform QC, merge summary statistics and perform meta-analysis.<br/>
Meta-analysis can be performed with `.meta()` with inverse-variance weighted or samplesize-weighted methods.<br/>
GWAMA as described in [Baselmans, et al. (2019)](https://www.nature.com/articles/s41588-018-0320-8) can be performed using the `.gwama()` function in merged summary statistics. <br/>
The plotting package uses matplotlib.pyplot for generating figures, so the functions are generally compatible with matplotlib.pyplot colors, and Figure and Axis objects. <br/>
Warning: merging with low_memory enabled is still highly experimental. <br/>

# Reference

Using the pysumstats package for a publication, or something similar? That is **awesome**! <br/>
There is no publication attached to this package, 
and I am not going to force anyone to reference me or make me a co-author or whatever, I want this to remain easily accessible. 
But I would greatly appreciate it if you add a link to this github, or a reference to it in the acknowledgements or something like that. <br/>
If you have any questions, want to help add methods or want to let me know you are planning a publication with this, you can get in touch via the [pypi website of this project](https://pypi.org/project/pysumstats/). <br/>
If you use the `.gwama()` method, please refer to the original publication: [Baselmans, et al. (2019)](https://www.nature.com/articles/s41588-018-0320-8).

# Installation

This package was made for Python 3.7. Clone the package directly from this github, or install with 

`pip3 install --upgrade pysumstats`


# Usage

`import pysumstats as sumstats`
###### Reading files
`s1 = sumstats.SumStats("sumstats1.csv.gz", phenotype='GWASsummary1')`
###### Reading data without sample size column: you will manually have to specify gwas sample size
`s2 = sumstats.SumStats("sumstats2.txt.gz", phenotype='GWASsummary2', gwas_n=350492)`
###### Reading data with column names not automatically recognized:
```
s3 = sumstats.SumStats("sumstats3.csv", phenotype='GWASsummary3',
                              column_names={
                                    'rsid': 'weird_name_for_rsid',
                                    'chr': 'weird_name_for_chr',
                                    'bp': 'weird_name_for_bp',
                                    'ea': 'weird_name_for_ea',
                                    'oa': 'weird_name_for_oa',
                                    'maf': 'weird_name_for_maf',
                                    'b': 'weird_name_for_b',
                                    'se': 'weird_name_for_se',
                                    'p': 'weird_name_for_p',
                                    'hwe': 'weird_name_for_p_hwe',
                                    'info': 'weird_name_for_info',
                                    'n': 'weird_name_for_n',
                                    'eaf': 'weird_name_for_eaf',
                                    'oaf': 'weird_name_for_oaf'})
```
###### Performing qc
```
s1.qc(maf=.01)
s2.qc(maf=.01, hwe=1e-6, info=.9)
s3.qc()  # MAF .01 is the default
```
###### Merging sumstats, low_memory option is still experimental so be careful with that
`merge1 = s1.merge(s2)`

###### Meta analysis
```
n_weighted_meta = merge1.meta_analyze(name='meta1', method='samplesize')  # N-weighted meta analysis
ivw_meta = merge1.meta_analyze(name='meta1', method='ivw')  # Standard inverse-variance weighted meta analysis
gwama = merge1.gwama(name='meta1', method='ivw')  # GWAMA as described in Baselmans, et al. (2019)
```
###### Additionally supports adding SNP heritabilities as weights
`exc_meta = exc.gwama(h2_snp={'ntr_exc': .01, 'ukb_ssoe': .02}, name='exc', method='ivw')`
###### And your own covariance matrix (called cov_Z in most R scripts)
```
# Either read it from a file:
import pandas as pd
cov_z = pd.read_csv('my_cov_z.csv') # Note it should be pandas dataframe with column names and index names equal to your phenotypes

# Or generate it from a phenotype file yourself:
phenotypes = pd.read_csv('my_phenotype_file.csv')
cov_z = sumstats.cov_matrix_from_phenotype_file(phenotypes, phenotypes=['GWASsummary1', 'GWASsummary2'])

gwama = exc.gwama(cov_matrix=cov_z, h2_snp={'GWASsummary1': .01, 'GWASsummary2': .02}, name='meta1', method='ivw')
```
###### See a summary of the result
`gwama.describe()`
###### See head of the data
`gwama.head()`
###### See head of all chromosomes
`gwama.head(n_chromosomes=23)`

###### QQ and Manhattan plots of the result
```
gwama.manhattan(filename='meta_manhattan.png')
gwama.qqplot(filename='meta_qq.png')
``` 

###### Save the result as csv
`exc.save('exc_sumstats.csv')`
###### Save the result as a pickle file (way faster to save and load back into Python)
`exc.save('exc_sumstats.pickle')`

###### Merge gwama results with another file:
`merged = gwama.merge(s3)`
###### Save prepped files for MR analysis in R:
```
merged.prep_for_mr(exposure='GWASsummary3', outcome='meta1',
                   filename=['GWAS3-Meta.csv', 'Meta-GWAS3.csv'],
                   p_cutoff=5e-8, bidirectional=True, index=False)
```
The resulting files will have the following column names, per specification of the MendelianRandomization package in R:

`rsid	chr	bp	exposure.A1	exposure.A2	outcome.A1	outcome.A2	exposure.se	exposure.b	outcome.se	outcome.b`

###### Some other stuff:
```
# See column names of the file
gpc_neuro.columns

# SumStats support for standard indexing is growing:
exc[0]  # Get the full output of the first SNP
exc[:10]  # Get the full output of the first 10 SNPs
exc[:10, 'p']  # Get the p value of the first 10 SNPs
exc['p']  # Get the p values of all SNPs
exc['rs78948828']  # Get the full output of 1 specific rsid
exc[['rs78948828', 'rs6057089', 'rs55957973']]  # Get the full output of multiple specific rsids
exc[['rs78948828', 'rs6057089', 'rs55957973'], 'p']  # Get the p-value for specific rsids

# If for whatever reason you want to do stuff with each SNP individually you can also loop over the entire file
for snp_output in exc:
    if exc['p'] < 5e-8:
        print('Yay significant SNP!')
    # do something


# If you only want to loop over some specific columns, you can
for rsid, b, se, p in exc[['rsid', 'b', 'se', 'p']].values:
    if p < 5e-8:
        print('Yay significant SNP!')


```

