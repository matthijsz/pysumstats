from scipy.stats import norm
import pandas as pd
import warnings
import os
import numpy as np
from .base import _BaseSumStats, _H5SSConnection



class MergedSumStats(_BaseSumStats):
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
                if not inplace:
                    self.data = merged_data
                    self.columns = self.data[1].columns
                    self.suffixes = ['x', 'y']
                    self._allign(ynames=['_{}'.format(x) for x in other.pheno_names])
                else:
                    newmergedss = MergedSumStats(data=merged_data, phenotypes=self.pheno_names + other.pheno_names,
                                                 merge_info={}, xy=['x', 'y'], variables=self.variables,
                                                 low_ram=self.low_ram, tmpdir=self.tmpdir, allign=False)
                    newmergedss._allign(ynames=['_{}'.format(x) for x in other.pheno_names])

    def describe(self, columns=None, per_chromosome=False):
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
