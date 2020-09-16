import pandas as pd
import os
import pickle
import numpy as np
from collections.abc import Iterable
from scipy.stats import norm
from .base import _BaseSumStats, _H5SSConnection
from .utils import _MultiWindowPlot, _recompute_maf
from .plot import (manhattan, qqplot, afplot, zzplot, pzplot)
from .exceptions import sumstatswarn

pd.options.mode.chained_assignment = None


class SumStats(_BaseSumStats):
    """Class for summary statistics of a single GWAS.

    :param path: Path to the file containing summary statistics. Should be a csv, or tab-delimited txt-file (.gz supported).
    :type path: str
    :param phenotype: Phenotype name
    :type phenotype: str
    :param gwas_n: Optional N subjects in the GWAS, for N-based meta analysis (if there is no N column in the summary statistics)
    :type gwas_n: int
    :param column_names: Optional dictionary of column names, if these are not automatically recognised. Keys should be: ['rsid', 'chr', 'bp', 'ea', 'oa', 'maf', 'b', 'se', 'p', 'hwe', 'info', 'n', 'eaf', 'oaf']
    :type column_names: dict
    :param data: Dataset for the new SumStats object, in general, don't specify this.
    :type data: dict
    :param low_ram: Whether to use the low_ram option for this SumStats object. Use this only when running into MemoryErrors. Enabling this option will read/write data from local storage rather then RAM. It will save lots of RAM, but it will significantly decrease processing speed.
    :type low_ram: bool
    :param tmpdir: Which directory to store the temporary files if low_ram is enabled.
    :type tmpdir: str
    :param kwargs: other keyword arguments to be passed to pandas.read_csv() method

    """

    def __init__(self, path, phenotype=None, gwas_n=None, column_names=None, data=None, low_ram=False,
                 tmpdir='sumstats_temporary', **kwargs):
        assert isinstance(path, str) or (path is None),  "path should be str"
        assert isinstance(phenotype, str) or (phenotype is None), "phenotype should be str"
        assert isinstance(gwas_n, int) or gwas_n is None, "gwas_n should be int"
        assert isinstance(low_ram, bool), "low_ram should be True or False"
        assert isinstance(tmpdir, str), "tmpdir should be str"
        super().__init__()
        self.low_ram = low_ram
        self.tmpdir = tmpdir
        self.plot_funcs = {'manhattan': self.manhattan, 'qqplot': self.qqplot, 'pzplot': self.pzplot}
        if phenotype is None:
            self.phenotype_name = path.split('/')[-1].split('.')[0]
        else:
            self.phenotype_name = phenotype
        self.column_names = column_names
        if path is not None:
            if path.endswith('.pickle'):
                with open(path, 'rb') as f:
                    sumstatsobj = pickle.load(f)
                if sumstatsobj.low_ram:
                    if not os.path.isfile(sumstatsobj.data.path):
                        raise ImportError(
                            'The data that was stored in {} using low_ram does not exist anymore.'.format(path))
                if self.__version__ != sumstatsobj.__version__:
                    if self.__version__.split('.')[0] != sumstatsobj.__version__.split('.')[0]:
                        pass  # What to do in case of major version differences: No problem yet.
                    else:
                        sumstatswarn('Imported sumstats version {} does not match current version {}'.format(
                            sumstatsobj.__version__, self.__version__))
                for k, v in {m: getattr(sumstatsobj, m) for m in dir(sumstatsobj) if
                             (not m.startswith('__'))}.items():
                    setattr(self, k, v)
                self.__version__ = sumstatsobj.__version__
            else:
                if 'sep' not in kwargs.keys():
                    if path.endswith('.txt') or path.endswith('.txt.gz') or path.endswith('.tsv') or path.endswith('.tsv.gz'):
                            self.data = pd.read_csv(path, sep='\t', **kwargs)
                    elif path.endswith('.csv') or path.endswith('.csv.gz'):
                        self.data = pd.read_csv(path, **kwargs)
                    else:
                        raise ImportError('Only .txt(.gz), .tsv(.gz) or .csv(.gz) files allowed when \'sep\' is not specified.')
                else:
                    self.data = pd.read_csv(path, **kwargs)
            self.gwas_n = gwas_n

            self._sync_columns()
            self._split()
        else:
            self.phenotype_name = phenotype
            self.data = data
            self.reset_index()
        self.qc_result = {}
        self.columns = self.data[1].columns

    def _sync_columns(self):
        """Internal function to synchronize column names and types with the defaults.

        :return: None

        """
        self.data.columns = [x.lower() for x in self.data.columns]
        if (len(self.data.columns)) != len(list(set(self.data.columns))):
            raise KeyError('Duplicate column names not allowed! (not case-sensitive)')
        required = ['rsid', 'chr', 'bp', 'ea', 'oa', 'maf', 'b', 'se', 'p']
        column_name_variants = {
            'rsid': ['rsid', 'rs', 'rsnumber'],
            'chr': ['chr', 'chromosome'],
            'bp': ['bp', 'pos'],
            'ea': ['a1', 'ea', 'effectallele', 'allele1', 'ref', 'refallele'],
            'oa': ['a2', 'oa', 'otherallele', 'allele0', 'allele2', 'alt', 'altallele'],
            'maf': ['maf', 'minorallelefrequency', 'minfreq'],
            'b': ['b', 'beta', 'effect'],
            'se': ['se', 'stderr'],
            'p': ['p', 'pval', 'pvalue', 'p_bolt_lmm_inf', 'p-value'],
            'hwe': ['p_hwe', 'hwe', 'hwe_p'],
            'info': ['info', 'rsq', 'r2'],
            'n': ['n', 'nchrobs', 'samplesize', 'totalsamplesize']
        }
        column_name_variants_af = {
            'eaf': ['a1freq', 'a1frq', 'a1f', 'eaf', 'freq1'],
            'oaf': ['a2freq', 'a2frq', 'a2f', 'oaf', 'freq2'],
        }
        column_types = {'rsid': str, 'chr': np.uint8, 'bp': np.uint32, 'ea': str, 'oa': str, 'maf': np.float32,
                        'b': np.float64, 'se': np.float64, 'p': np.float64, 'hwe': np.float32, 'info': np.float32,
                        'n': np.uint32, 'eaf': np.float32, 'oaf': np.float32}
        if self.column_names is not None:
            for k, v in self.column_names:
                if (k not in column_name_variants) and (k not in column_name_variants_af):
                    raise KeyError('{} not recognized as column type'.format(k))
                elif k in column_name_variants:
                    column_name_variants[k] = [v]
                elif k in column_name_variants_af:
                    column_name_variants_af[k] = [v]
        found, found_af, not_found = [], [], []
        for target, options in column_name_variants.items():
            for opt in options:
                if opt in self.data.columns:
                    self.data.rename(columns={opt: target}, inplace=True)
                    found += [target]
                    break
            if (target == 'rsid') and ('rsid' not in found):
                for opt in ['snpname', 'snp', 'snpid', 'cptid', 'markername']:
                    if opt in self.data.columns:
                        self.data.rename(columns={opt: target}, inplace=True)
                        found.append(target)
                        sumstatswarn(
                            'Assuming column {} represents rsid as no other rsid column was found.'.format(opt))
                        break
        if self.data['rsid'].isnull().values.any():
            sumstatswarn('Missing RSIDs were replaced with CHR:BP')
            midx = self.data['rsid'].isnull()
            self.data.loc[midx, 'rsid'] = self.data.loc[midx, 'chr'].map(str) + ':' + self.data.loc[midx, 'bp'].map(str)
        for target, options in column_name_variants_af.items():
            for opt in options:
                if opt in self.data.columns:
                    self.data.rename(columns={opt: target}, inplace=True)
                    found_af.append(target)
        if 'maf' not in found:
            if len(found_af) == 0:
                raise KeyError('No allele frequency column found')
            if ('eaf' in found_af) and ('oaf' not in found_af):
                self.data['oaf'] = 1 - self.data['eaf']
            elif ('eaf' not in found_af) and ('oaf' in found_af):
                self.data['eaf'] = 1 - self.data['oaf']
            self.data['maf'] = self.data[['eaf', 'oaf']].min(axis=1)
            found.append('maf')
        for target in required:
            if target not in found:
                not_found.append(target)
        if 'n' not in self.data.columns:
            if self.gwas_n is None:
                sumstatswarn('No N found or specified, samplesize based meta analysis not possible')
            else:
                self.data['n'] = self.gwas_n
        else:
            if np.any(self.data['n'].isna()):
                nnna = np.sum(self.data['n'].isna())
                sumstatswarn('{} missing values in sample size column were mean imputed'.format(nnna))
                self.data.loc[self.data['n'].isna(), 'n'] = self.data['n'].mean()
        if np.any(self.data['bp'].isna()):
            naidx = np.where(self.data['bp'].isna())
            self.data.loc[naidx[0], 'bp'] = (self.data.loc[naidx[0] - 1, 'bp'] + 1).tolist()
            if np.any(self.data['bp'].isna()):
                raise ValueError('Too many missing values in basepair column')
            else:
                sumstatswarn('{} missing values in basepair column were imputed to be (bp at the previous row)+1'.format(len(naidx)))
        if np.any(self.data['chr'].isna()):
            naidx = np.where(self.data['chr'].isna())
            self.data.loc[naidx[0], 'chr'] = (self.data.loc[naidx[0] - 1, 'chr']).tolist()
            if np.any(self.data['chr'].isna()):
                raise ValueError('Too many missing values in chromosome column')
            else:
                sumstatswarn('{} missing values in chromosome column were imputed to be the same as chromsome of the previous row'.format(len(naidx)))
        try:
            self.data['chr'] = self.data['chr'].astype(int)
        except ValueError:
            self.data['chr'] = self.data['chr'].str.replace('X', '23')
            try:
                self.data['chr'] = self.data['chr'].astype(int)
            except ValueError:
                raise ValueError('Could not convert chromosome column to integer')
        if len(not_found) > 0:
            raise KeyError('Could not find columns: {}'.format(', '.join(not_found)))
        self.data['ea'] = self.data['ea'].str.upper()
        self.data['oa'] = self.data['oa'].str.upper()
        self.variables = list(self.data.columns)
        for column, targettype in column_types.items():
            if column in self.data.columns:
                self.data[column].astype(targettype, copy=False)

    def _split(self):
        """ Internal function to split the initial large dataset into chunks based on chromosome

        :return: None

        """
        if not self.low_ram:
            new_data = {}
        else:
            new_data = _H5SSConnection(self.phenotype_name, tmpdir=self.tmpdir)
        for c in range(1, 24):
            new_data[c] = self.data.loc[self.data['chr'] == c, :].copy()
            if len(new_data[c]) == 0:
                sumstatswarn('No data found for chromosome {}'.format(c))
        self.data = new_data

    def qc(self, maf=None, hwe=None, info=None, **kwargs):
        """Basic GWAS quality control function.

        :param maf: Minor allele frequency cutoff, will drop SNPs where MAF < cutoff. Default: 0.1
        :type maf: float or None
        :param hwe: Hardy-Weinberg Equilibrium cutoff, will drop SNPs where HWE < cutoff, if specified and HWE column is present in the data.
        :type hwe: float or None
        :param info: Imputation quality cutoff, will drop SNPs where Info < cutoff, if specified and Info column is present in the data.
        :type info: float or None
        :param kwargs: Other columns to filter on, keyword should be column name, SNPs whill be dropped where the value < argument.
        :return: None

        """
        assert isinstance(maf, float) or (maf is None), "maf should be float"
        assert isinstance(hwe, float) or (hwe is None), "hwe should be float"
        assert isinstance(info, float) or (info is None), "info should be float"
        qc_vals = dict(maf=.01)
        qc_info = dict(org_len=0, maf=0, new_len=0)
        if maf is not None:
            qc_vals['maf'] = maf
        if hwe is not None:
            if 'hwe' not in self.variables:
                sumstatswarn('HWE qc cutoff specified but HWE column was not found, skipping QC step.')
            else:
                qc_vals['hwe'] = hwe
                qc_info['hwe'] = 0
        if info is not None:
            if 'info' not in self.variables:
                sumstatswarn('Info qc cutoff specified but info column was not found, skipping QC step.')
            else:
                qc_vals['info'] = info
                qc_info['info'] = 0
        for k, v in kwargs.items():
            if k not in self.variables:
                raise KeyError('Column {} not found.'.format(k))
            else:
                qc_vals[k] = v
                qc_info[k] = 0
        for c in range(1, 24):
            data = self.data[c]
            qc_info['org_len'] += len(data)
            for var, cutoff in qc_vals.items():
                qc_info[var] += len(data.loc[data[var] >= cutoff, :])
            for var, cutoff in qc_vals.items():
                data = data.loc[data[var] >= cutoff, :].copy()
            qc_info['new_len'] += len(data)
            self.data[c] = data
        self.qc_result = qc_info

    def merge(self, other, how='inner', low_memory=False):
        """Merge with other SumStats object(s).

        :param other: Other sumstats object, or list of other SumStats objects.
        :type other: :class:`pysumstats.plot.SumStats` or list
        :param how: Type of merge.
        :type how: str
        :param low_memory: Enable to use a more RAM-efficient merging method (WARNING: still untested)
        :type low_memory: bool
        :return: :class:`pysumstats.plot.MergedSumStats` object

        """
        assert isinstance(other, SumStats) or isinstance(other, list), "other should be SumStats or list"
        if isinstance(other, list):
            for o in other:
                assert isinstance(o, SumStats), "items in other should be SumStats"
            merged = self.merge(other[0], low_memory=low_memory)
            merged.merge(other[1:], inplace=True, low_memory=low_memory)
            return merged
        else:
            if self.phenotype_name == other.phenotype_name:
                sumstatswarn('Phenotype names were equal, converting to {0}_x and {0}_y'.format(self.phenotype_name))
                self.phenotype_name = self.phenotype_name + '_x'
                other.phenotype_name = other.phenotype_name + '_y'
            merge_info = dict(xlen=0, ylen=0, newlen=0)
            joined = [x for x in self.variables if x in other.variables]
            for x in [x for x in self.variables if x not in joined]:
                sumstatswarn(
                    'Could not find column {} in pysumstats for {}, column dropped'.format(x, other.phenotype_name))
            for x in [x for x in other.variables if x not in joined]:
                sumstatswarn(
                    'Could not find column {} in pysumstats for {}, column dropped'.format(x, self.phenotype_name))
            if not self.low_ram:
                merged_data = {}
            else:
                merged_data = _H5SSConnection(
                    '{}{}{}{}{}'.format(self.phenotype_name[0], self.phenotype_name[-1], other.phenotype_name[0],
                                        other.phenotype_name[1], len(os.listdir('sumstats_temporary'))),
                    tmpdir=self.tmpdir)
            if 'n' not in self.columns:
                sumstatswarn('No sample size column found for {}, using 1 for calculation of overall MAF'.format(self.phenotype_name))
            if 'n' not in other.columns:
                sumstatswarn('No sample size column found for {}, using 1 for calculation of overall MAF'.format(other.phenotype_name))
            if low_memory and (not how == 'inner'):
                raise NotImplementedError('low_memory only allows for an inner merge.')
            for c in range(1, 24):
                data_x = self.data[c].copy()
                merge_info['xlen'] += len(data_x)
                data_y = other.data[c].copy()
                merge_info['ylen'] += len(data_y)
                data_m_x = data_x.rename(columns={k: '{}_{}'.format(k, self.phenotype_name) for k in data_x.columns if
                                                  k not in ['rsid', 'chr', 'bp']})
                data_m_y = data_y.rename(columns={k: '{}_{}'.format(k, other.phenotype_name) for k in data_y.columns if
                                                  k not in ['rsid', 'chr', 'bp']})
                if low_memory:
                    if np.any(data_m_x.duplicated(subset='rsid')) or np.any(data_m_y.duplicated(subset='rsid')):
                        raise KeyError('Duplicated rsids in either dataset will cause low_memory merge to fail')
                    merged_data[c] = data_m_x.loc[data_m_x['rsid'].isin(data_m_y['rsid'].tolist()), :]
                    data_m_y = data_m_y.loc[data_m_y['rsid'].isin(data_m_x['rsid'].tolist()), :]
                    merged_data[c].sort_values(by='rsid', inplace=True)
                    data_m_y.sort_values(by='rsid', inplace=True)
                    merged_data[c].reset_index(inplace=True, drop=True)
                    data_m_y.reset_index(inplace=True, drop=True)
                    data_m_y = data_m_y[
                        [x + '_' + other.phenotype_name for x in joined if x not in ['rsid', 'chr', 'bp']]]
                    for x in data_m_y.columns:
                        merged_data[c][x] = data_m_y[x]
                        data_m_y.drop(axis=1, labels=[x], inplace=True)
                    merged_data[c].sort_values(by='bp', inplace=True)
                else:
                    joined_x = ['rsid', 'chr', 'bp'] + [x + '_' + self.phenotype_name for x in joined if
                                                        x not in ['rsid', 'chr', 'bp']]
                    joined_y = ['rsid'] + [x + '_' + other.phenotype_name for x in joined if
                                           x not in ['rsid', 'chr', 'bp']]
                    merged_data[c] = data_m_x[joined_x].merge(data_m_y[joined_y], on='rsid', how=how)
                merged_data[c]['maf'] = _recompute_maf(merged_data[c], [self.phenotype_name, other.phenotype_name])
                merge_info['newlen'] += len(merged_data[c])
            return MergedSumStats(data=merged_data, phenotypes=[self.phenotype_name, other.phenotype_name],
                                  merge_info=merge_info, variables=joined,
                                  xy=[self.phenotype_name, other.phenotype_name], low_ram=self.low_ram,
                                  tmpdir=self.tmpdir)

    def describe(self, columns=None, per_chromosome=False):
        """Get a summary of the data.

        :param columns: List of column names to print summary for (default: ['b', 'se', 'p'])
        :type columns: list.
        :param per_chromosome: Enable to return a list of summary dataframes per chromosome
        :type per_chromosome: bool.
        :return: pd.Dataframe, or list

        """
        assert isinstance(columns, list) or (columns is None), "columns should be a list"
        assert isinstance(per_chromosome, bool), "per_chromosome should be True or False"
        if columns is None:
            columns = ['b', 'se', 'p']
        if (not isinstance(columns, list)) and isinstance(columns, str):
            columns = [columns]
        sum_cols = []
        for c in columns:
            if c in self.data[1].columns:
                sum_cols.append(c)
            else:
                sumstatswarn('{} not found in columns, skipping.'.format(c))
        return self._get_summary(sum_cols, per_chromosome)

    def manhattan(self, **kwargs):
        """Generate a manhattan plot using this sumstats data

        :param kwargs: keyworded arguments to be passed to :func:`pysumstats.plot.manhattan`
        :return: None, or (fig, ax)

        """
        manhattan(self[['rsid', 'chr', 'bp', 'p']], **kwargs)

    def qqplot(self, **kwargs):
        """Generate a QQ-plot using this sumstats data

        :param kwargs: keyworded arguments to be passed to :func:`pysumstats.plot.qqplot`
        :return: None, or (fig, ax)

        """
        qqplot(self['p'].values, **kwargs)

    def pzplot(self, **kwargs):
        """Generate a PZ-plot using this sumstats data

        :param kwargs: keyworded arguments to be passed to :func:`pysumstats.plot.pzplot`
        :return: None, or (fig, ax)

        """
        pzplot(self[['b', 'se', 'p']], **kwargs)


class MergedSumStats(_BaseSumStats):
    """Class containing merged summary statistics. In general you will not create a MergedSumStats object manually.

    :param data: dataset containing merged summary statistics
    :type data: dict
    :param phenotypes: list of phenotype names.
    :type phenotypes: list
    :param merge_info: Dict with information on the merge
    :type merge_info: dict
    :param variables: list of variables contained in the data.
    :type variables: list
    :param xy: x and y suffixes (to be used in _allign)
    :type xy: list
    :param low_ram: Whether to use the low_ram option for this MergedSumStats object (passed down from SumStats). Use this only when running into MemoryErrors. Enabling this option will read/write data from local storage rather then RAM. It will save lots of RAM, but it will gratly decrease processing speed.
    :type low_ram: bool
    :param tmpdir: Which directory to store the temporary files if low_ram is enabled (passed down from SumStats).
    :type tmpdir: str
    :param allign: Enable to auto-allign SNPs
    :type allign: bool

    """

    def __init__(self, data, phenotypes, merge_info, variables, xy, low_ram=False, tmpdir='sumstats_temporary',
                 allign=True):
        super().__init__()
        self.pheno_names = phenotypes
        self.data = data
        self.info = merge_info
        self.variables = variables
        self.suffixes = xy
        self.columns = self.data[1].columns
        if allign:
            self._allign()
        self.reset_index()
        self.phenotype_name = None
        self.low_ram = low_ram
        self.tmpdir = tmpdir
        self.plot_funcs = {'manhattan': self.manhattan, 'qqplot': self.qqplot, 'pzplot': self.pzplot,
                           'afplot': self.afplot, 'zzplot': self.zzplot}

    def _allign(self, ynames=None):
        """Function to allign SNPs to the first phenotype.

        :param ynames: Optional argument of multiple phenotypes that should be alligned.
        :type ynames: list
        :return: None

        """
        if self.suffixes[0] is None:
            suffix_x = ''
        else:
            suffix_x = '_' + self.suffixes[0]
        if self.suffixes[1] is None:
            suffix_y = ''
        else:
            suffix_y = '_' + self.suffixes[1]
        eax, oax = 'ea{}'.format(suffix_x), 'oa{}'.format(suffix_x)
        if ynames is None:
            eay, oay, bys = 'ea{}'.format(suffix_y), 'oa{}'.format(suffix_y), ['b{}'.format(suffix_y)]
            eafys, oafys = ['eaf{}'.format(suffix_y)], ['oaf{}'.format(suffix_y)]
        else:
            eay, oay, bys = 'ea{}'.format(suffix_y), 'oa{}'.format(suffix_y), ['b{}'.format(y) for y in ynames]
            eafys, oafys = ['eaf{}'.format(y) for y in ynames], ['oaf{}'.format(y) for y in ynames]
        dropped = 0
        for c in self.data.keys():
            data = self.data[c].copy()
            org_len = len(data)
            for a in [eax, oax, eay, oay]:
                data[a] = data[a].astype(str)
            flip = data.loc[(data[eax] == data[oay]) & (data[oax] == data[eay]), bys[0]].index
            for by in bys:
                data.loc[flip, by] *= -1
            data.loc[flip, eay] = data.loc[flip, eax]
            data.loc[flip, oay] = data.loc[flip, oax]
            for fy in eafys + oafys:
                if fy in data.columns:
                    data.loc[flip, fy] = abs(1 - data.loc[flip, fy])
            data2 = data.loc[data[eax] == data[eay], :].copy()
            data2 = data2.loc[data2[oax] == data2[oay], :].copy()
            dropped += (org_len - len(data))
            data2.drop(axis=1, labels=[eay, oay], inplace=True)
            data2.rename(columns={eax: 'ea', oax: 'oa'}, inplace=True)
            self.data[c] = data2
        self.columns = self.data[1].columns
        if dropped > 0:
            sumstatswarn('Dropped {} SNPs due to allele mismatch'.format(dropped))
        self.suffixes = [None, None]

    def meta_analyze(self, name='meta', method='ivw'):
        """Meta analyze all GWAS summary statistics contained in this object.

        :param name: New phenotype name to use for the new SumStats object (default: 'meta')
        :type name: str
        :param method: Meta-analysis method to use, should be one of ['ivw', 'samplesize'], default: 'ivw'
        :type method: str
        :return: :class:`pysumstats.SumStats` object.

        """
        assert isinstance(name, str), "name should be str"
        if not self.low_ram:
            new_data = {}
        else:
            new_data = _H5SSConnection(name, tmpdir=self.tmpdir)
        if method not in ['ivw', 'samplesize']:
            raise KeyError('method should be one of [\'ivw\', \'samplesize\']')
        columns = self.data[1].columns
        missing_n = []
        for p in self.pheno_names:
            if 'n_{}'.format(p) not in columns:
                missing_n.append(p)
        if len(missing_n) > 0:
            if method == 'samplesize':
                raise KeyError('Missing sample size column for {}.'.format(', '.join(missing_n)))
            else:
                sumstatswarn('Missing sample size column for {}, output N will be incorrect.'.format(', '.join(missing_n)))
        for c in self.data.keys():
            data = self.data[c]
            n_dat = data[['rsid', 'chr', 'bp', 'ea', 'oa', 'maf']]
            if method == 'samplesize':
                for p in self.pheno_names:
                    data.loc[:, 'z_{}'.format(p)] = data['b_{}'.format(p)] / data['se_{}'.format(p)]
                    data.loc[:, 'w_{}'.format(p)] = np.sqrt(data['n_{}'.format(p)])
                    data.loc[:, 'z_{}'.format(p)] = data['z_{}'.format(p)] * data['w_'.format(p)]
                zsums = data[['z_'.format(x) for x in self.pheno_names]].sum(axis=1)
                wsums = data[['w_'.format(x) for x in self.pheno_names]].sum(axis=1)
                n_dat.loc[:, 'z'] = zsums / np.sqrt(np.sum(wsums ** 2))
            else:
                for p in self.pheno_names:
                    data.loc[:, 'w_{}'.format(p)] = 1 / (data.loc[:, 'se_{}'.format(p)] ** 2)
                    data.loc[:, 'bw_{}'.format(p)] = data.loc[:, 'b_{}'.format(p)] * data.loc[:, 'w_{}'.format(p)]
                bwsums = data[['bw_{}'.format(x) for x in self.pheno_names]].sum(axis=1)
                wsums = data[['w_{}'.format(x) for x in self.pheno_names]].sum(axis=1)
                n_dat.loc[:, 'se'] = np.sqrt(1 / wsums)
                n_dat.loc[:, 'b'] = bwsums / wsums
                n_dat.loc[:, 'z'] = n_dat.loc[:, 'b'] / n_dat.loc[:, 'se']
            n_dat.loc[:, 'p'] = norm.sf(abs(n_dat.loc[:, 'z'])) * 2
            n_dat['n'] = 0
            for p in self.pheno_names:
                if 'n_{}'.format(p) in n_dat.columns:
                    n_dat['n'] += data['n_{}'.format(p)]
            new_data[c] = n_dat
        new_columns = list(new_data[c].columns)
        new_data = SumStats(path=None, phenotype=name, data=new_data)
        new_data.variables = new_columns
        return new_data

    def gwama(self, cov_matrix=None, h2_snp=None, name='gwama'):
        """Multivariate meta analysis as described in Baselmans, et al. 2019.

        :param cov_matrix: Covariance matrix, defaults to generating a correlation matrix of Z-scores
        :type cov_matrix: pd.Dataframe
        :param h2_snp: Dict of SNP heritabilities per GWAS, to use as additional weights. Defaults to all 1's.
        :type h2_snp: dict
        :param name: New phenotype name to use in the new SumStats object (default: 'gwama')
        :return: :class:`pysumstats.SumStats` object

        """
        assert isinstance(cov_matrix, pd.DataFrame), "cov_matrix should be a pd.DataFrame"
        assert isinstance(name, str), "name should be a string"
        if h2_snp is None:
            sumstatswarn('h2-snp not specified, using ones instead. This will bias the estimates')
            h2_snp = {x: 1 for x in self.pheno_names}
        elif (not isinstance(h2_snp, pd.DataFrame)) or isinstance(h2_snp, dict) or (not isinstance(h2_snp, pd.Series)):
            raise KeyError('h2-snp should be a dataframe or dictionary with phenotype names as keys or column names.')
        else:
            for p in self.pheno_names:
                try:
                    _ = h2_snp[p]
                except KeyError:
                    raise KeyError('Could not find h2_snp value for {}'.format(p))
        if cov_matrix is None:
            cov_matrix = pd.DataFrame(index=self.pheno_names, columns=self.pheno_names)
            vars = {}
            for p in self.pheno_names:
                vars[p] = []
                cov_matrix.loc[p, p] = 1
                for c in range(1, 24):
                    data = self.data[c]
                    data['z_{}'.format(p)] = data['b_{}'.format(p)] / data['se_{}'.format(p)]
                    vars[p].append(self.data[1]['z_{}'.format(p)].var() * (len(data) / len(self)))
                    self.data[c] = data
            for p1, p2 in [(i, i2) for n, i in enumerate(self.pheno_names) for i2 in self.pheno_names[(n + 1):]]:
                tmp_statistics = {}
                for c in range(1, 24):
                    if len(self.data[c]) > 0:
                        covar = self.data[c][['z_{}'.format(p1), 'z_{}'.format(p2)]].cov().iloc[0, 1] * (
                                len(self.data[c]) / len(self))
                        tmp_statistics[c] = dict(p1_var=vars[p1][c], p2_var=vars[p2][c], covar=covar)
                tmp_df = pd.DataFrame.from_dict(tmp_statistics)
                cov_matrix.loc[p1, p2] = tmp_df.loc['covar', :].sum() / (
                        np.sqrt(tmp_df.loc['p1_var', :].sum()) * np.sqrt(tmp_df.loc['p2_var', :].sum()))
                cov_matrix.loc[p2, p1] = cov_matrix.loc[p1, p2]
        else:
            if not isinstance(cov_matrix, pd.DataFrame):
                raise KeyError(
                    'Incorrect cov_matrix specified, should be pandas dataframe with phenotypes as columns and indices.')
            if (not sorted(list(cov_matrix.index)) == sorted(list(self.pheno_names))) or (
                    not sorted(list(cov_matrix.columns)) == sorted(list(self.pheno_names))):
                raise IndexError(
                    'Incorrect cov_matrix specified, should be pandas dataframe with phenotypes as columns and indices.')
        if not self.low_ram:
            new_data = {}
        else:
            new_data = _H5SSConnection(name, tmpdir=self.tmpdir)
        i_cov_matrix = pd.DataFrame(np.linalg.inv(cov_matrix.to_numpy(dtype=np.float64)), index=cov_matrix.index,
                                    columns=cov_matrix.columns)
        for c in range(1, 24):
            n_dat = self.data[c][['rsid', 'chr', 'bp', 'ea', 'oa', 'maf']].copy()
            data = self.data[c]
            for p in self.pheno_names:
                data['z_{}'.format(p)] = data['b_{}'.format(p)] / data['se_{}'.format(p)]
                data['w_{}'.format(p)] = np.sqrt(data['n_{}'.format(p)] * h2_snp[p])
                data['wz_{}'.format(p)] = data['w_{}'.format(p)] * data['z_{}'.format(p)]
            data['cw'] = np.array(
                [data['w_{}'.format(p)] * data['w_{}'.format(p2)] * cov_matrix.loc[p, p2]
                 for p in self.pheno_names for p2 in [x for x in self.pheno_names]]).sum(axis=0)
            n_dat['z'] = data[['wz_{}'.format(x) for x in self.pheno_names]].sum(axis=1) / np.sqrt(
                data['cw'])
            n_dat.loc[:, 'p'] = norm.sf(abs(n_dat.loc[:, 'z'])) * 2
            n_dat['n'] = np.array([np.sqrt(data['n_{}'.format(p)]) * np.sqrt(
                data['n_{}'.format(p2)]) * i_cov_matrix.loc[p, p2] for p in self.pheno_names for p2 in
                                   [x for x in self.pheno_names]]).sum(axis=0)
            n_obs = data[['n_{}'.format(p) for p in self.pheno_names]].sum(axis=1)
            n_dat['eaf'] = np.array(
                [data['n_{}'.format(p)] * data['eaf_{}'.format(p)] for p in self.pheno_names]).sum(
                axis=0) / n_obs
            n_dat['maf'] = pd.concat([n_dat['eaf'], (1 - n_dat['eaf'])]).min(axis=0)
            n_dat['b'] = n_dat['z'] / np.sqrt(n_obs * 2 * n_dat['eaf'] * (1 - n_dat['eaf']))
            n_dat['se'] = (1 / np.sqrt(n_obs)) * (1 / np.sqrt(2 * n_dat['eaf'] * (1 - n_dat['eaf'])))
            new_data[c] = n_dat
        new_columns = list(new_data[1].columns)
        new_data = SumStats(path=None, phenotype=name, data=new_data)
        new_data.variables = new_columns
        return new_data

    def merge(self, other, inplace=False, how='inner', low_memory=False):
        """ Merge with other SumStats or MergedSumstats object(s).

        :param other: :class:`pysumstats.SumStats`, or :class:`pysumstats.MergedSumStats` object, or a list of SumStats, MergedSumstats objects
        :type other: :class:`pysumstats.SumStats`, :class:`pysumstats.MergedSumStats`, or list.
        :param inplace: Enable to store the new data in the current MergedSumStats object. (currently not supported when low_ram is enabled)
        :type inplace: bool
        :param how: Type of merge, for now only implemented for merges with :class:`pysumstats.SumStats` objects
        :type how: str
        :param low_memory: Enable to use a more RAM-efficient merging method (WARNING: still untested)
        :type low_memory: bool
        :return: None, or :class:`pysumstats.MergedSumStats` object.

        """
        assert isinstance(other, SumStats) or isinstance(other, MergedSumStats) or isinstance(other, list), "other should be SumStats, MergedSumStats or list"
        assert isinstance(inplace, bool), "inplace should be True or False"
        assert isinstance(how, str), "how should be a string"
        assert isinstance(low_memory, bool), "low_memory should be True or False"
        if isinstance(other, list):
            for o in other:
                assert isinstance(o, SumStats) or isinstance(o, MergedSumStats), "items in other should be SumStats or MergedSumStats"
        if self.low_ram and inplace:
            sumstatswarn('Inplace merging with low_ram currently does not yield any performance benefits.')
        if isinstance(other, list):
            if inplace:
                for o in other:
                    self.merge(o, True, low_memory)
            else:
                newmerged = self.merge(other[0], False, low_memory)
                for o in other[1:]:
                    newmerged.merge(o, True, low_memory)
                return newmerged
        else:
            if low_memory and (not how == 'inner'):
                raise NotImplementedError('low_memory only allows for an inner merge.')
            if other.phenotype_name is not None:
                if other.phenotype_name in self.pheno_names:
                    sumstatswarn(
                        'Phenotype name already in the dataset, converting to and {0}2'.format(other.phenotype_name))
                    other.phenotype_name = other.phenotype_name + '2'
                self.suffixes[1] = other.phenotype_name
                merge_info = dict(xlen=0, ylen=0, newlen=0)
                joined = [x for x in other.variables if x in self.variables]
                for x in [x for x in other.variables if x not in joined]:
                    sumstatswarn(
                        'Could not find column {} in pysumstats for {}, column dropped'.format(x, other.phenotype_name))
                merged_data = {}
                for c in range(1, 24):
                    data_x = self.data[c].copy()
                    merge_info['xlen'] += len(data_x)
                    data_y = other.data[c].copy()
                    merge_info['ylen'] += len(data_y)
                    data_m_y = data_y.rename(
                        columns={k: '{}_{}'.format(k, other.phenotype_name) for k in data_y.columns if
                                 k not in ['rsid', 'chr', 'bp']})
                    if low_memory:
                        if np.any(data_x.duplicated(subset='rsid')) or np.any(data_m_y.duplicated(subset='rsid')):
                            raise KeyError('Missing rsids in either dataset will cause low_memory merge to fail')
                        merged_data[c] = data_x.loc[data_x['rsid'].isin(data_m_y['rsid'].tolist()), :]
                        data_m_y = data_m_y.loc[data_m_y['rsid'].isin(data_x['rsid'].tolist()), :]
                        del data_y
                        merged_data[c].sort_values(by='rsid', inplace=True)
                        data_m_y.sort_values(by='rsid', inplace=True)
                        merged_data[c].reset_index(inplace=True, drop=True)
                        data_m_y.reset_index(inplace=True, drop=True)
                        data_m_y = data_m_y[
                            [x + '_' + other.phenotype_name for x in joined if x not in ['rsid', 'chr', 'bp']]]
                        for x in data_m_y.columns:
                            merged_data[c][x] = data_m_y[x]
                            data_m_y.drop(axis=1, labels=[x], inplace=True)
                        merged_data[c].sort_values(by='bp', inplace=True)
                    else:
                        joined_y = ['rsid'] + [x + '_' + other.phenotype_name for x in joined if
                                               x not in ['rsid', 'chr', 'bp']]
                        merged_data[c] = data_x.merge(data_m_y[joined_y], on='rsid', how=how)
                    merged_data[c]['maf'] = _recompute_maf(merged_data[c], self.pheno_names + [other.phenotype_name])
                    merge_info['newlen'] += len(merged_data[c])
                    new_phenos = self.pheno_names + [other.phenotype_name]
                if not inplace:
                    return MergedSumStats(data=merged_data, phenotypes=new_phenos, merge_info=merge_info,
                                          variables=joined, xy=self.suffixes, low_ram=self.low_ram, tmpdir=self.tmpdir)
                else:
                    self.__init__(merged_data, new_phenos, merge_info, joined, self.suffixes, self.low_ram)
            else:
                for n, other_name in enumerate(other.pheno_names):
                    if other_name in self.pheno_names:
                        other.pheno_names[n] = other_name + '_y'
                        for c in range(1, 24):
                            for x in other.variables:
                                odat = other.data[c]
                                odat.rename(columns={x + '_' + other_name: x + '_' + other_name + '_y'}, inplace=True)
                                other.data[c] = odat
                if not self.low_ram:
                    merged_data = {}
                else:
                    tmpfilename = ''.join([''.join([j[0], j[1]]) for j in self.pheno_names]) + ''.join(
                        [''.join([j[0], j[1]]) for j in other.pheno_names])
                    merged_data = _H5SSConnection(
                        '{}{}'.format(tmpfilename, len(os.listdir('sumstats_temporary'))), tmpdir=self.tmpdir)
                for c in range(1, 24):
                    sdat = self.data[c]
                    odat = other.data[c]
                    s_id = pd.DataFrame(sdat['rsid'])
                    o_id = pd.DataFrame(odat['rsid'])
                    s_id['s_idx'] = s_id.index
                    o_id['o_idx'] = o_id.index
                    idm = s_id.merge(o_id, on='rsid')
                    sdat.rename(columns={'ea': 'ea_x', 'oa': 'oa_x'}, inplace=True)
                    odat.rename(columns={'ea': 'ea_y', 'oa': 'oa_y'}, inplace=True)
                    other_cols = [x for x in odat.columns if x not in ['rsid', 'chr', 'bp']]
                    sdat = sdat.loc[idm['s_idx'], :].reset_index(drop=True)
                    odat = odat.loc[idm['o_idx'], other_cols].reset_index(drop=True)
                    for oc in odat.columns:
                        sdat[oc] = odat[oc]
                    merged_data[c] = sdat
                    merged_data[c]['maf'] = _recompute_maf(merged_data[c], self.pheno_names + other.pheno_names)
                if inplace:
                    self.data = merged_data
                    self.columns = self.data[1].columns
                    self.suffixes = ['x', 'y']
                    self._allign(ynames=['_{}'.format(x) for x in other.pheno_names])
                else:
                    newmergedss = MergedSumStats(data=merged_data, phenotypes=self.pheno_names + other.pheno_names,
                                                 merge_info={}, xy=['x', 'y'], variables=self.variables,
                                                 low_ram=self.low_ram, tmpdir=self.tmpdir, allign=False)
                    newmergedss._allign(ynames=['_{}'.format(x) for x in other.pheno_names])
                    return newmergedss

    def describe(self, columns=None, per_chromosome=False):
        """ Get a summary of the data.

        :param columns: List of column names to print summary for (default: ['b', 'se', 'p'])
        :type columns: list
        :param per_chromosome: Enable to return a list of summary dataframes per chromosome
        :type per_chromosome: bool
        :return: pd.Dataframe, or list of pd.Dataframes

        """
        assert isinstance(columns, list) or (columns is None), "columns should be a list"
        assert isinstance(per_chromosome, bool), "per_chromosome should be True or False"
        if columns is None:
            columns = ['b', 'se', 'p']
        if (not isinstance(columns, list)) and isinstance(columns, str):
            columns = [columns]
        sum_cols = []
        for c in columns:
            if c not in self.data[1].columns:
                if '{}_{}'.format(c, self.pheno_names[0]) in self.data[1].columns:
                    for c2 in ['{}_{}'.format(c, p) for p in self.pheno_names]:
                        sum_cols.append(c2)
                else:
                    sumstatswarn('column {} not found'.format(c))
            else:
                sum_cols.append(c)
        return self._get_summary(sum_cols, per_chromosome)

    def prep_for_mr(self, exposure, outcome, filename=None, p_cutoff=None, bidirectional=False, **kwargs):
        """Save a pre-formatted file to use with the MendelianRandomization package in R.

        :param exposure: phenotype name to use as exposure.
        :type exposure: str
        :param outcome: phenotype name to use as outcome.
        :type outcome: str
        :param filename: Path to where the resulting file(s) should be stored, or list of paths if bidirectional=True
        :type filename: str, list or None
        :param p_cutoff: Optional p-value cut-off to apply. Will include SNPs where P > p_cutoff
        :type p_cutoff: float
        :param bidirectional: Enable to store two files (exposure=exposure, outcome=outcome), and (exposure=outcome, outcome=exposure)
        :type bidirectional: bool
        :param kwargs: Additional keyword arguments to be passed to pandas to_csv function.
        :return: None

        """
        assert isinstance(exposure, str), "exposure should be a string"
        assert isinstance(outcome, str), "exposure should be a string"
        assert isinstance(filename, str) or isinstance(filename, list) or (filename is None), "filename should be a string or list of strings"
        assert isinstance(p_cutoff, float), "p_cutoff should be a float"
        assert isinstance(bidirectional, bool), "bidirectional should be True or False"
        if isinstance(filename, list):
            for o in filename:
                assert isinstance(o, str), "filename should be a string or list of strings"
        if bidirectional and (not isinstance(filename, list)):
            raise ValueError(
                'If bidirectional, filename should be a list of filenames: (exp-outcome) name and (outcome-exp) name.')
        mr_data = []
        if exposure not in self.pheno_names:
            raise KeyError('{} not found in phenotypes'.format(exposure))
        if outcome not in self.pheno_names:
            raise KeyError('{} not found in phenotypes'.format(outcome))
        for c in self.data.keys():
            if p_cutoff is None:
                mr_data_ = self.data[c].copy()
            else:
                mr_data_ = self.data[c].loc[self.data[c]['p_{}'.format(exposure)] <= p_cutoff, :].copy()
            newcols = [['se', 'se'], ['b', 'b']]
            newcols2 = [[exposure, 'exposure'], [outcome, 'outcome']]
            renames = {'{}_{}'.format(old, x): '{}.{}'.format(x2, new) for x, x2 in newcols2 for old, new in newcols}
            mr_data_.rename(columns=renames, inplace=True)
            mr_data_.rename(columns={'ea': 'exposure.A1', 'oa': 'exposure.A2'}, inplace=True)
            mr_data_['outcome.A1'] = mr_data_['exposure.A1']
            mr_data_['outcome.A2'] = mr_data_['exposure.A2']
            mr_data.append(mr_data_[
                               ['rsid', 'chr', 'bp', 'exposure.A1', 'exposure.A2', 'outcome.A1', 'outcome.A2'] + list(
                                   renames.values())])
        mr_data = pd.concat(mr_data, axis=0)
        if filename is None:
            mr_data.to_csv('MR-prepped_{}_{}.csv'.format(exposure, outcome), **kwargs)
        else:
            if not bidirectional:
                mr_data.to_csv(filename, **kwargs)
            else:
                mr_data.to_csv(filename[0], **kwargs)
        if bidirectional:
            if filename is None:
                self.prep_for_mr(exposure=outcome, outcome=exposure, filename=None, p_cutoff=p_cutoff,
                                 bidirectional=False, **kwargs)
            else:
                self.prep_for_mr(exposure=outcome, outcome=exposure, filename=filename[1], p_cutoff=p_cutoff,
                                 bidirectional=False, **kwargs)

    def manhattan(self, filename=None, phenotypes=None, nrows=None, ncols=None, figsize=None, dpi=300, **kwargs):
        """Generates manhattan plots for each phenotype (or specified phenotypes) in merged GWAS data.

        :param filename: Target file to save the resulting figure to (if no name is specified, fig and axes are returned)
        :type filename: str
        :param phenotypes: List of phenotype names to plot manhattans for (defaults to plotting all phenotypes)
        :type phenotypes: list
        :param nrows: Specify number of rows in the figure ( defaults to int(ceil(len(phenotypes)/ncols)) )
        :type nrows: int
        :param ncols: Specify number of columns in the figure ( defaults to int(log2(len(phenotypes)/2)) )
        :type ncols: int
        :param figsize: Specify width and height of figure in inches ( defaults to (ncols*8, nrows*4) )
        :type figsize: (int, int)
        :param dpi: DPI setting to use when saving the figure.
        :type dpi: int
        :param kwargs: Other keyword arguments to be passed to :func:`pysumstats.plot.manhattan`
        :return: (fig, axes) or None.
        """
        assert isinstance(filename, str) or (filename is None), "filename should be a string or None"
        assert isinstance(phenotypes, list) or (phenotypes is None), "phenotypes should be a list of strings or None"
        assert isinstance(dpi, int), "dpi should be an integer"
        if phenotypes is None:
            phenotypes = self.pheno_names
        for p in phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        plotwindow = _MultiWindowPlot(len(phenotypes), nrows=nrows, ncols=ncols, figsize=figsize, shape='rect')
        dat = self[['rsid', 'chr', 'bp']]
        for phenotype in phenotypes:
            fig, ax = plotwindow.get_next_ax()
            dat['p'] = self['p_{}'.format(phenotype)]
            if 'title' not in kwargs.keys():
                manhattan(dat, fig=fig, ax=ax, title=phenotype, **kwargs)
            else:
                manhattan(dat, fig=fig, ax=ax, **kwargs)
        return plotwindow.finish(filename=filename, dpi=dpi)

    def qqplot(self, filename=None, phenotypes=None, nrows=None, ncols=None, figsize=None, dpi=300, **kwargs):
        """Generates QQ-plots for each phenotype (or specified phenotypes) in merged GWAS data.

        :param filename: Target file to save the resulting figure to (if no name is specified, fig and axes are returned)
        :type filename: str
        :param phenotypes: List of phenotype names to plot QQ-plots for (defaults to plotting all phenotypes)
        :type phenotypes: list
        :param nrows: Specify number of rows in the figure ( defaults to int(sqrt(len(phenotypes))) )
        :type nrows: int
        :param ncols: Specify number of columns in the figure ( defaults to int(ceil(len(phenotypes)/nrows)) )
        :type ncols: int
        :param figsize: Specify width and height of figure in inches ( defaults to (ncols*5, nrows*5) )
        :type figsize: (int, int)
        :param dpi: DPI setting to use when saving the figure.
        :type dpi: int
        :param kwargs: Other keyword arguments to be passed to :func:`pysumstats.plot.qqplot`
        :return: (fig, axes) or None.

        """
        assert isinstance(filename, str) or (filename is None), "filename should be a string or None"
        assert isinstance(phenotypes, list) or (phenotypes is None), "phenotypes should be a list of strings or None"
        assert isinstance(dpi, int), "dpi should be an integer"
        if phenotypes is None:
            phenotypes = self.pheno_names
        for p in phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        plotwindow = _MultiWindowPlot(len(phenotypes), nrows=nrows, ncols=ncols, figsize=figsize, shape='square')
        for phenotype in phenotypes:
            fig, ax = plotwindow.get_next_ax()
            if 'title' not in kwargs.keys():
                qqplot(self['p_{}'.format(phenotype)].values, fig=fig, ax=ax, title=phenotype, **kwargs)
            else:
                qqplot(self['p_{}'.format(phenotype)].values, fig=fig, ax=ax, **kwargs)
        return plotwindow.finish(filename=filename, dpi=dpi)

    def pzplot(self, filename=None, phenotypes=None, nrows=None, ncols=None, figsize=None, dpi=300, **kwargs):
        """Generates PZ-plots for each phenotype (or specified phenotypes) in merged GWAS data.

        :param filename: Target file to save the resulting figure to (if no name is specified, fig and axes are returned)
        :type filename: str
        :param phenotypes: List of phenotype names to plot PZ-plots for (defaults to plotting all phenotypes)
        :type phenotypes: list
        :param nrows: Specify number of rows in the figure ( defaults to int(sqrt(len(phenotypes))) )
        :type nrows: int
        :param ncols: Specify number of columns in the figure ( defaults to int(ceil(len(phenotypes)/nrows)) )
        :type ncols: int
        :param figsize: Specify width and height of figure in inches ( defaults to (ncols*5, nrows*5) )
        :type figsize: (int, int)
        :param dpi: DPI setting to use when saving the figure.
        :type dpi: int
        :param kwargs: Other keyword arguments to be passed to :func:`pysumstats.plot.pzplot`
        :return: (fig, axes) or None.

        """
        assert isinstance(filename, str) or (filename is None), "filename should be a string or None"
        assert isinstance(phenotypes, list) or (phenotypes is None), "phenotypes should be a list of strings or None"
        assert isinstance(dpi, int), "dpi should be an integer"
        if phenotypes is None:
            phenotypes = self.pheno_names
        for p in phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        plotwindow = _MultiWindowPlot(len(phenotypes), nrows=nrows, ncols=ncols, figsize=figsize, shape='square')
        for phenotype in phenotypes:
            dat = self[['{}_{}'.format(v, phenotype) for v in ['b', 'se', 'p']]].copy()
            dat.rename(columns={'{}_{}'.format(v, phenotype): v for v in ['b', 'se', 'p']}, inplace=True)
            fig, ax = plotwindow.get_next_ax()
            if 'title' not in kwargs.keys():
                pzplot(dat, fig=fig, ax=ax, title='PZ_{}'.format(phenotype), **kwargs)
            else:
                pzplot(dat, fig=fig, ax=ax, **kwargs)
        return plotwindow.finish(filename=filename, dpi=dpi)

    def afplot(self, ref_phenotypes=None, other_phenotypes=None, filename=None, nrows=None, ncols=None, figsize=None,
               dpi=300, **kwargs):
        """Generates AF comparison plots for merged GWAS data.

        :param ref_phenotypes: List of phenotypes to use as reference (defaults to all phenotypes)
        :type ref_phenotypes: list
        :param other_phenotypes: List of phenotypes to compare reference to (defaults to all phenotypes, overlapping plots will be dropped)
        :type other_phenotypes: list
        :param filename: Target file to save the resulting figure to (if no name is specified, fig and axes are returned)
        :type filename: str
        :param nrows: Specify number of rows in the figure ( defaults to int(ceil(n_plots/ncols)) )
        :type nrows: int
        :param ncols: Specify number of columns in the figure ( defaults to int(sqrt(n_plots)))
        :type ncols: int
        :param figsize: Specify width and height of figure in inches ( defaults to (ncols*5, nrows*5) )
        :type figsize: (int, int)
        :param dpi: DPI setting to use when saving the figure.
        :type dpi: int
        :param kwargs: Other keyword arguments to be passed to :func:`pysumstats.plot.manhattan`
        :return: (fig, axes) or None.

        """
        assert isinstance(filename, str) or (filename is None), "filename should be a string or None"
        assert isinstance(ref_phenotypes, list) or (ref_phenotypes is None), "ref_phenotypes should be a string or None"
        assert isinstance(other_phenotypes, list) or (other_phenotypes is None), "other_phenotypes should be a string or None"
        assert isinstance(dpi, int), "dpi should be an integer"
        if ref_phenotypes is None:
            ref_phenotypes = self.pheno_names
        for p in ref_phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        if other_phenotypes is None:
            other_phenotypes = self.pheno_names
        for p in other_phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        n_plots = len(ref_phenotypes) * len(other_phenotypes)
        plotwindow = _MultiWindowPlot(n_plots, nrows=nrows, ncols=ncols, figsize=figsize, shape='square')
        for p1 in ref_phenotypes:
            eafp1 = self['eaf_{}'.format(p1)].values
            for p2 in other_phenotypes:
                fig, ax = plotwindow.get_next_ax()
                if p1 != p2:
                    afplot(eafp1, self['eaf_{}'.format(p2)].values,
                           refname=r'$EAF_{}$'.format('{' + p1.replace('_', r'\_') + '}'),
                           othername=r'$EAF_{}$'.format('{' + p2.replace('_', r'\_') + '}'), fig=fig, ax=ax,
                           **kwargs)
        return plotwindow.finish(filename=filename, dpi=dpi)

    def zzplot(self, ref_phenotypes=None, other_phenotypes=None, filename=None, nrows=None, ncols=None, figsize=None,
               dpi=300, **kwargs):
        """Generates ZZ comparison plots for merged GWAS data.

        :param ref_phenotypes: List of phenotypes to use as reference (defaults to all phenotypes)
        :type ref_phenotypes: list
        :param other_phenotypes: List of phenotypes to compare reference to (defaults to all phenotypes, overlapping plots will be dropped)
        :type other_phenotypes: list
        :param filename: Target file to save the resulting figure to (if no name is specified, fig and axes are returned)
        :type filename: str
        :param nrows: Specify number of rows in the figure ( defaults to int(ceil(n_plots/ncols)) )
        :type nrows: int
        :param ncols: Specify number of columns in the figure ( defaults to int(sqrt(n_plots)))
        :type ncols: int
        :param figsize: Specify width and height of figure in inches ( defaults to (ncols*5, nrows*5) )
        :type figsize: (int, int)
        :param dpi: DPI setting to use when saving the figure.
        :type dpi: int
        :param kwargs: Other keyword arguments to be passed to :func:`pysumstats.plot.zzplot`
        :return: (fig, axes) or None.

        """
        assert isinstance(filename, str) or (filename is None), "filename should be a string or None"
        assert isinstance(ref_phenotypes, list) or (ref_phenotypes is None), "ref_phenotypes should be a string or None"
        assert isinstance(other_phenotypes, list) or (other_phenotypes is None), "other_phenotypes should be a string or None"
        assert isinstance(dpi, int), "dpi should be an integer"
        if ref_phenotypes is None:
            ref_phenotypes = self.pheno_names
        for p in ref_phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        if other_phenotypes is None:
            other_phenotypes = self.pheno_names
        for p in other_phenotypes:
            assert p in self.pheno_names, "{} not found in phenotypes".format(p)
        n_plots = len(ref_phenotypes) * len(other_phenotypes)
        plotwindow = _MultiWindowPlot(n_plots, nrows=nrows, ncols=ncols, figsize=figsize, shape='square')
        for p1 in ref_phenotypes:
            zz1 = self[['b_{}'.format(p1), 'se_{}'.format(p1)]].copy()
            zz1.rename(columns={'b_{}'.format(p1): 'b', 'se_{}'.format(p1): 'se'}, inplace=True)
            for p2 in other_phenotypes:
                fig, ax = plotwindow.get_next_ax()
                if p1 != p2:
                    zz2 = self[['b_{}'.format(p2), 'se_{}'.format(p2)]].copy()
                    zz2.rename(columns={'b_{}'.format(p2): 'b', 'se_{}'.format(p2): 'se'}, inplace=True)
                    zzplot(zz1, zz2, xname=r'$Z_{}$'.format('{' + p1.replace('_', r'\_') + '}'),
                           yname='r$Z_{}$'.format('{' + p2.replace('_', r'\_') + '}'),
                           fig=fig,
                           ax=ax, **kwargs)
        return plotwindow.finish(filename=filename, dpi=dpi)
