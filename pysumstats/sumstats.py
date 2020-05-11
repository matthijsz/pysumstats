import pandas as pd
import warnings
import os
import numpy as np
from scipy.stats import norm
from .base import _BaseSumStats, _H5SSConnection
from .plot import manhattan
from .plot import qqplot

pd.options.mode.chained_assignment = None


class SumStats(_BaseSumStats):
    """Class for summary statistics of a single GWAS.

    :param path: Path to the file containing summary statistics. Should be a csv, or tab-delimited txt-file (.gz supported).
    :type path: str.
    :param phenotype: Phenotype name
    :type phenotype: str.
    :param gwas_n: Optional N subjects in the GWAS, for N-based meta analysis (if there is no N column in the summary statistics)
    :type gwas_n: int
    :param column_names: Optional dictionary of column names, if these are not automatically recognised. Keys should be: ['rsid', 'chr', 'bp', 'ea', 'oa', 'maf', 'b', 'se', 'p', 'hwe', 'info', 'n', 'eaf', 'oaf']
    :type column_names: dict.
    :param data: Dataset for the new SumStats object, in general, don't specify this.
    :type data: dict.
    :param low_ram: Whether to use the low_ram option for this SumStats object. Use this only when running into MemoryErrors. Enabling this option will read/write data from local storage rather then RAM. It will save lots of RAM, but it will gratly decrease processing speed.
    :type low_ram: bool.
    :param tmpdir: Which directory to store the temporary files if low_ram is enabled.
    :type tmpdir: str.

    """
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
        """Internal function to synchronize column names with the defaults.

        :return: None. Renamed data is stored inplace.

        """
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
        """ Internal function to split the initial large dataset into chunks based on chromosome

        :return: None. Split data is stored inplace.

        """
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
        """Basic GWAS quality control function.

        :param maf: Minor allele frequency cutoff, will drop SNPs where MAF < cutoff. Default: 0.1
        :param hwe: Hardy-Weinberg Equilibrium cutoff, will drop SNPs where HWE < cutoff, if specified and HWE column is present in the data.
        :param info: Imputation quality cutoff, will drop SNPs where Info < cutoff, if specified and Info column is present in the data.
        :return: None. Data is stored inplace

        """
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
        """Merge with other SumStats object(s).

        :param other: Other sumstats object, or list of other SumStats objects.
        :type other: SumStats or list.
        :param low_memory: Enable to use a more RAM-efficient merging method (WARNING: still untested)
        :type low_memory: bool.
        :return: MergedSumstats object.

        """
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
        """Get a summary of the data.

        :param columns: List of column names to print summary for (default: ['b', 'se', 'p'])
        :type columns: list.
        :param per_chromosome: Enable to return a list of summary dataframes per chromosome
        :type per_chromosome: bool.
        :return: pd.Dataframe, or list

        """
        if columns is None:
            columns = ['b', 'se', 'p']
        if (not isinstance(columns, list)) and isinstance(columns, str):
            columns = [columns]
        sum_cols = []
        for c in columns:
            if c in self.data[1].columns:
                sum_cols.append(c)
        return self._get_summary(sum_cols, per_chromosome)

    def manhattan(self, **kwargs):
        """Generate a manhattan plot using this sumstats data

        :param kwargs: keyworded arguments to be passed to :func:`pysumstats.plot.manhattan`
        :return: None, or (fig, ax)

        """
        manhattan(self[['rsid', 'chr', 'bp', 'p']], **kwargs)


class MergedSumStats(_BaseSumStats):
    """Class containing merged summary statistics. In general you will not create a MergedSumStats object manually.

    :param data: dataset containing merged summary statistics
    :type data: dict.
    :param phenotypes: list of phenotype names.
    :type phenotypes: list.
    :param merge_info: Dict with information on the merge
    :type merge_info: dict.
    :param variables: list of variables contained in the data.
    :type variables: list
    :param xy: x and y suffixes (to be used in _allign)
    :type xy: list
    :param low_ram: Whether to use the low_ram option for this MergedSumStats object (passed down from SumStats). Use this only when running into MemoryErrors. Enabling this option will read/write data from local storage rather then RAM. It will save lots of RAM, but it will gratly decrease processing speed.
    :type low_ram: bool
    :param tmpdir: Which directory to store the temporary files if low_ram is enabled (passed down from SumStats).
    :type tmpdir: str.
    :param allign: Enable to auto-allign SNPs
    :type allign: bool.

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

    def _allign(self, ynames=None):
        """Function to allign SNPs to the first phenotype.

        :param ynames: Optional argument of multiple phenotypes that should be alligned.
        :type ynames: list.
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
            warnings.warn('Dropped {} SNPs due to allele mismatch'.format(dropped))
        self.suffixes = [None, None]

    def meta_analyze(self, name='meta', method='ivw'):
        """Meta analyze all GWAS summary statistics contained in this object.

        :param name: New phenotype name to use for the new SumStats object (default: 'meta')
        :type name: str.
        :param method: Meta-analysis method to use, should be one of ['ivw', 'samplesize'], default: 'ivw'
        :type method: str.
        :return: :class:`pysumstats.SumStats` object.

        """
        if not self.low_ram:
            new_data = {}
        else:
            new_data = _H5SSConnection(name, tmpdir=self.tmpdir)
        if method not in ['ivw', 'samplesize']:
            raise KeyError('method should be one of [\'ivw\', \'samplesize\']')
        for c in self.data.keys():
            data = self.data[c]
            n_dat = data[['rsid', 'chr', 'bp', 'ea', 'oa']]
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
            new_data[c] = n_dat
        new_columns = list(new_data[c].columns)
        new_data = SumStats(path=None, phenotype=name, data=new_data)
        new_data.variables = new_columns
        return new_data

    def gwama(self, cov_matrix=None, h2_snp=None, name='gwama'):
        """Multivariate meta analysis as described in Baselmans, et al. 2019.

        :param cov_matrix: Covariance matrix, defaults to generating a correlation matrix of Z-scores
        :type cov_matrix: pd.Dataframe.
        :param h2_snp: Dict of SNP heritabilities per GWAS, to use as additional weights. Defaults to all 1's.
        :type h2_snp: dict.
        :param name: New phenotype name to use in the new SumStats object (default: 'gwama')
        :return: :class:`pysumstats.SumStats` object.

        """
        if h2_snp is None:
            warnings.warn('h2-snp not specified, using ones instead. This will bias the estimates')
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
            n_dat = self.data[c][['rsid', 'chr', 'bp', 'ea', 'oa']].copy()
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

    def merge(self, other, inplace=False, low_memory=False):
        """ Merge with other SumStats or MergedSumstats object(s).

        :param other: :class:`pysumstats.SumStats`, or :class:`pysumstats.MergedSumStats` object, or a list of SumStats, MergedSumstats objects
        :type other: :class:`pysumstats.SumStats`, :class:`pysumstats.MergedSumStats`, or list.
        :param inplace: Enable to store the new data in the current MergedSumStats object. (currently not supported when low_ram is enabled)
        :type inplace: bool.
        :param low_memory: Enable to use a more RAM-efficient merging method (WARNING: still untested)
        :type low_memory: bool.
        :return: :class:`pysumstats.MergedSumStats` object.

        """
        if self.low_ram and inplace:
            return NotImplementedError('Merging inplace when low_ram is enabled is currently not supported.')
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
            if other.phenotype_name is not None:
                if other.phenotype_name in self.pheno_names:
                    warnings.warn(
                        'Phenotype name already in the dataset, converting to and {0}2'.format(other.phenotype_name))
                    other.phenotype_name = other.phenotype_name + '2'
                self.suffixes[1] = other.phenotype_name
                merge_info = dict(xlen=0, ylen=0, newlen=0)
                joined = [x for x in other.variables if x in self.variables]
                for x in [x for x in other.variables if x not in joined]:
                    warnings.warn(
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
                        merged_data[c] = data_x.merge(data_m_y[joined_y], on='rsid')
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
        :type columns: list.
        :param per_chromosome: Enable to return a list of summary dataframes per chromosome
        :type per_chromosome: bool.
        :return: pd.Dataframe, or list of pd.Dataframes containing data summary

        """
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
                sum_cols.append(c)
        return self._get_summary(sum_cols, per_chromosome)

    def prep_for_mr(self, exposure, outcome, filename=None, p_cutoff=None, bidirectional=False, **kwargs):
        """Save a pre-formatted file to use with the MendelianRandomization package in R.

        :param exposure: phenotype name to use as exposure.
        :type exposure: str.
        :param outcome: phenotype name to use as outcome.
        :type outcome: str.
        :param filename: Path to where the resulting file(s) should be stored, or list of paths is bidirectional=True
        :type filename: str or list.
        :param p_cutoff: Optional p-value cut-off to apply. Will include SNPs where P > p_cutoff
        :type p_cutoff: float.
        :param bidirectional: Enable to store two files (exposure=exposure, outcome=outcome), and (exposure=outcome, outcome=exposure)
        :type bidirectional: bool.
        :param kwargs: Additional keyword arguments to be passed to pandas to_csv function.
        :return: None, resulting data is stored on local storage.

        """
        if bidirectional and (not isinstance(filename, list)):
            raise ValueError(
                'If bidirectional, filename should be a list of filenames: (exp-outcome)-name and (outcome-exp) name.')
        mr_data = []
        for c in self.data.keys():
            if p_cutoff is None:
                mr_data_ = self.data[c].copy()
            else:
                mr_data_ = self.data[c].loc[self.data[c]['p_{}'.format(exposure)] <= p_cutoff, :]
            newcols = [['se', 'se'], ['b', 'b'], ['se', 'se']]
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
