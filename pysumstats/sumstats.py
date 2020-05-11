import pandas as pd
import warnings
import os
import numpy as np
from .base import _BaseSumStats, _H5SSConnection
from .mergedsumstats import MergedSumStats

pd.options.mode.chained_assignment = None


class SumStats(_BaseSumStats):
    '''
    Class for summary statistics of a single GWAS.
    :param path: Path to the file containing summary statistics. Should be a csv, or tab-delimited txt-file (.gz supported)
    :type path: str.
    :param phenotype: Phenotype name
    :type phenotype: str.
    :param gwas_n: Optional N subjects in the GWAS, for N-based meta analysis (if there is no N column in the summary statistics)
    :type gwas_n: int
    :param column_names: Optional dictionary of column names, if these are not automatically recognised. Keys should be:
    ['rsid', 'chr', 'bp', 'ea', 'oa', 'maf', 'b', 'se', 'p', 'hwe', 'info', 'n', 'eaf', 'oaf']
    :type column_names: dict.
    :param data: Dataset for the new SumStats object, in general, don't specify this.
    :type data: dict.
    :param low_ram: Whether to use the low_ram option for this SumStats object. Use this only when running into MemoryErrors.
    Enabling this option will read/write data from local storage rather then RAM. It will save lots of RAM, but it will gratly decrease processing speed.
    :type low_ram: bool.
    :param tmpdir: Which directory to store the temporary files if low_ram is enabled.
    :type tmpdir: str.
    '''
    def __init__(self, path, phenotype=None, gwas_n=None, column_names=None, data=None, low_ram=False,
                 tmpdir='sumstats_temporary'):
        super().__init__()
        self.low_ram = low_ram
        self.tmpdir = tmpdir
        if phenotype is None:
            self.phenotype_name = path.split('/')[-1].split('.')[0]
        else:
            self.phenotype_name = phenotype
        self.column_names = column_names
        if path is not None:
            if path.endswith('.txt') or path.endswith('.txt.gz'):
                self.data = pd.read_csv(path, sep='\t')
            elif path.endswith('.csv') or path.endswith('.csv.gz'):
                self.data = pd.read_csv(path)
            else:
                raise ImportError('Only .txt(.gz) and .csv(.gz) files are supported.')
            self.gwas_n = gwas_n
            self.qc_result = {}
            self._sync_colnames()
            self._split()
        else:
            self.phenotype_name = phenotype
            self.data = data
            self.reset_index()
        self.columns = self.data[1].columns

    def _sync_colnames(self):
        '''
        Internal function to synchronize column names with the defaults
        :return: None. Renamed data is stored inplace.
        '''
        self.data.columns = [x.lower() for x in self.data.columns]
        if (len(self.data.columns)) != len(list(set(self.data.columns))):
            raise KeyError('Duplicate column names not allowed! (not case-sensitive)')
        required = ['rsid', 'chr', 'bp', 'ea', 'oa', 'maf', 'b', 'se', 'p']
        column_name_variants = {
            'rsid': ['rsid', 'rs', 'rsnumber'],
            'chr': ['chr', 'chromosome'],
            'bp': ['bp', 'pos'],
            'ea': ['a1', 'ea', 'effectallele', 'allele1'],
            'oa': ['a2', 'oa', 'otherallele', 'allele0', 'allele2'],
            'maf': ['maf', 'minorallelefrequency', 'minfreq'],
            'b': ['b', 'beta', 'effect'],
            'se': ['se', 'stderr'],
            'p': ['p', 'pval', 'pvalue', 'p_bolt_lmm_inf', 'p-value'],
            'hwe': ['p_hwe', 'hwe', 'hwe_p'],
            'info': ['info', 'rsq', 'r2'],
            'n': ['n', 'nchrobs', 'samplesize', 'totalsamplesize']
        }
        column_name_variants2 = {
            'eaf': ['a1freq', 'a1frq', 'a1f', 'eaf', 'freq1'],
            'oaf': ['a2freq', 'a2frq', 'a2f', 'oaf', 'freq2'],
        }
        if self.column_names is not None:
            for k, v in self.column_names:
                if (k not in column_name_variants) and (k not in column_name_variants2):
                    raise KeyError('{} not recognized as column type'.format(k))
                elif k in column_name_variants:
                    column_name_variants[k] = [v]
                elif k in column_name_variants2:
                    column_name_variants2[k] = [v]
        found, found2, not_found = [], [], []
        for target, options in column_name_variants.items():
            for opt in options:
                if opt in self.data.columns:
                    self.data.rename(columns={opt: target}, inplace=True)
                    found += [target]
            if (target == 'rsid') and ('rsid' not in found):
                for opt in ['snp', 'snpid', 'cptid', 'markername']:
                    if opt in self.data.columns:
                        self.data.rename(columns={opt: target}, inplace=True)
                        found.append(target)
                        warnings.warn(
                            'Assuming column {} represents rsid as no other rsid column was found'.format(opt))
        if self.data['rsid'].isnull().values.any():
            warnings.warn('Missing RSIDs were replace with CHR:BP')
            self.data.loc[self.data['rsid'].isnull(), 'rsid'] = self.data.loc[self.data['rsid'].isnull(), 'chr'].map(
                str) + ':' + self.data.loc[self.data['rsid'].isnull(), 'bp'].map(str)
        for target, options in column_name_variants2.items():
            for opt in options:
                if opt in self.data.columns:
                    self.data.rename(columns={opt: target}, inplace=True)
                    found2.append(target)
        if 'maf' not in found:
            if len(found2) == 0:
                raise KeyError('No allele frequency column found')
            if ('eaf' in found2) and ('oaf' not in found2):
                self.data['oaf'] = 1 - self.data['eaf']
            elif ('eaf' not in found2) and ('oaf' in found2):
                self.data['eaf'] = 1 - self.data['oaf']
            self.data['maf'] = self.data[['eaf', 'oaf']].min(axis=1)
            found.append('maf')
        for target in required:
            if target not in found:
                not_found.append(target)
        if 'n' not in self.data.columns:
            if self.gwas_n is None:
                warnings.warn('No N found or specified, samplesize based meta analysis not possible')
            else:
                self.data['n'] = self.gwas_n
        if len(not_found) > 0:
            raise KeyError('Could not find columns: {}'.format(', '.join(not_found)))
        self.data['ea'] = self.data['ea'].str.upper()
        self.data['oa'] = self.data['oa'].str.upper()
        self.variables = list(self.data.columns)

    def _split(self):
        '''
        Internal function to split the initial large dataset into chunks based on chromosome
        :return: None. Split data is stored inplace.
        '''
        try:
            self.data['chr'] = self.data['chr'].astype(int)
        except ValueError:
            self.data['chr'].str.replace('X', '23')
            try:
                self.data['chr'] = self.data['chr'].astype(int)
            except ValueError:
                raise ValueError('Could not convert chromosome columnt to integer')
        if not self.low_ram:
            new_data = {}
        else:
            new_data = _H5SSConnection(self.phenotype_name, tmpdir=self.tmpdir)
        for c in range(1, 24):
            new_data[c] = self.data.loc[self.data['chr'] == c, :].copy()
            if len(new_data[c]) == 0:
                warnings.warn('No data found for chromosome {}'.format(c))
        self.data = new_data

    def qc(self, maf=None, hwe=None, info=None):
        '''
        Basic GWAS quality control function.
        :param maf: Minor allele frequency cutoff, will drop SNPs where MAF < cutoff. Default: 0.1
        :param hwe: Hardy-Weinberg Equilibrium cutoff, will drop SNPs where HWE < cutoff, if specified and HWE column is present in the data.
        :param info: Imputation quality cutoff, will drop SNPs where Info < cutoff, if specified and Info column is present in the data.
        :return: None. Data is stored inplace
        '''
        qc_vals = dict(maf=.01)
        qc_info = dict(org_len=0, maf=0, new_len=0)
        if maf is not None:
            qc_vals['maf'] = maf
        if hwe is not None:
            if 'hwe' not in self.variables:
                warnings.warn('HWE qc cutoff specified but HWE column was not found, skipping QC step.')
            else:
                qc_vals['hwe'] = hwe
                qc_info['hwe'] = 0
        if info is not None:
            if 'info' not in self.variables:
                warnings.warn('Info qc cutoff specified but info column was not found, skipping QC step.')
            else:
                qc_vals['info'] = info
                qc_info['info'] = 0
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

    def merge(self, other, low_memory=False):
        '''
        Merge with other SumStats object(s).
        :param other: Other sumstats object, or list of other SumStats objects.
        :type other: SumStats or list.
        :param low_memory: Enable to use a more RAM-efficient merging method (WARNING: still untested)
        :type low_memory: bool.
        :return: MergedSumstats object.
        '''
        if isinstance(other, list):
            merged = self.merge(other[0], low_memory=low_memory)
            merged.merge(other[1:], inplace=True, low_memory=low_memory)
            return merged
        else:
            if self.phenotype_name == other.phenotype_name:
                warnings.warn('Phenotype names were equal, converting to {0}_x and {0}_y'.format(self.phenotype_name))
                self.phenotype_name = self.phenotype_name + '_x'
                other.phenotype_name = other.phenotype_name + '_y'
            merge_info = dict(xlen=0, ylen=0, newlen=0)
            joined = [x for x in self.variables if x in other.variables]
            for x in [x for x in self.variables if x not in joined]:
                warnings.warn(
                    'Could not find column {} in pysumstats for {}, column dropped'.format(x, other.phenotype_name))
            for x in [x for x in other.variables if x not in joined]:
                warnings.warn(
                    'Could not find column {} in pysumstats for {}, column dropped'.format(x, self.phenotype_name))
            if not self.low_ram:
                merged_data = {}
            else:
                merged_data = _H5SSConnection(
                    '{}{}{}{}{}'.format(self.phenotype_name[0], self.phenotype_name[-1], other.phenotype_name[0],
                                        other.phenotype_name[1], len(os.listdir('sumstats_temporary'))),
                    tmpdir=self.tmpdir)
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
                    merged_data[c] = data_m_x[joined_x].merge(data_m_y[joined_y], on='rsid')
                merge_info['newlen'] += len(merged_data[c])
            return MergedSumStats(data=merged_data, phenotypes=[self.phenotype_name, other.phenotype_name],
                                  merge_info=merge_info, variables=joined,
                                  xy=[self.phenotype_name, other.phenotype_name], low_ram=self.low_ram,
                                  tmpdir=self.tmpdir)

    def describe(self, columns=None, per_chromosome=False):
        '''
        Get a summary of the data.
        :param columns: List of column names to print summary for (default: ['b', 'se', 'p'])
        :type columns: list.
        :param per_chromosome: Enable to return a list of summary dataframes per chromosome
        :type per_chromosome: bool.
        :return: pd.Dataframe, or list of pd.Dataframes containing data summary
        '''
        if columns is None:
            columns = ['b', 'se', 'p']
        if (not isinstance(columns, list)) and isinstance(columns, str):
            columns = [columns]
        sum_cols = []
        for c in columns:
            if c in self.data[1].columns:
                sum_cols.append(c)
        return self._get_summary(sum_cols, per_chromosome)
