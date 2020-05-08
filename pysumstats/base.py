import pandas as pd
import warnings
import pickle
import os
import numpy as np


class _Loc:
    def __init__(self, master):
        self.master = master

    def __getitem__(self, item):
        res = []
        for c in self.master.data.keys():
            res.append(self.master.data[c].loc[item])
        return pd.concat(res)

    def __setitem__(self, key, value):
        for c in self.master.data.keys():
            self.master.data[c].loc[key] = value


class _H5SSConnection:
    def __init__(self, filename, tmpdir):
        if not os.path.isdir(tmpdir):
            os.makedirs(tmpdir)
        self.path = '{}/{}'.format(tmpdir, filename)
        if os.path.isfile(self.path):
            suffix = 1
            while os.path.isfile(self.path + str(suffix)):
                suffix += 1
            self.path += str(suffix)


    def __getitem__(self, item):
        return pd.read_hdf(self.path, key='chr' + str(item))

    def __setitem__(self, key, value):
        value.to_hdf(self.path, key='chr' + str(key))

    def keys(self):
        return range(1, 24)

    def values(self):
        for i in self.keys():
            yield self[i]

    def items(self):
        for i in self.keys():
            yield i, self[i]

    def close(self):
        os.rename(self.path)


class _BaseSumStats:
    def __init__(self):
        self.loc = _Loc(self)

    def __len__(self):
        x = 0
        for data in self.data.values():
            x += len(data)
        return x

    def __getitem__(self, item):
        if isinstance(item, int):
            for c in list(range(1, 24)):
                if item > len(self.data[c]):
                    item -= len(self.data[c])
                else:
                    break
            if (c == 23) and (item >= len(self.data[23])):
                raise ValueError('Index not found')
            return self.data[c].loc[item, :]
        elif isinstance(item, slice):
            start, stop, step = item.start, item.stop, item.step
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)
            if step is None:
                step = 1
            return pd.concat([self[i] for i in range(start, stop, step)], axis=1).transpose()
        elif isinstance(item, list):
            if isinstance(item[0], str):
                if item[0] not in self.data[1].columns:
                    return pd.concat([self[i] for i in item], axis=0)
                else:
                    return pd.concat([self[i] for i in item], axis=1)
            else:
                return pd.concat([self[i] for i in item], axis=1).transpose()
        elif isinstance(item, str):
            if item in list(self.data[1].columns):
                return pd.concat([self.data[x][item] for x in range(1, 24)])
            else:
                for c in range(1, 24):
                    if item in self.data[c]['rsid'].tolist():
                        return self.data[c].loc[self.data[c]['rsid'] == item, :]
            raise ValueError('Index not found')
        elif isinstance(item, tuple):
            return self[item[0]][item[1]]
        else:
            raise ValueError('Index not supported')

    def __setitem__(self, key, value):
        if isinstance(value, int) or isinstance(value, float) or isinstance(value, str):
            for c in range(1, 24):
                data = self.data[c]
                data[key] = value
                self.data[c] = data
        if isinstance(value, pd.Series):
            if len(value) == len(self):
                start = 0
                for c in range(1, 24):
                    data = self.data[c]
                    data.loc[:, key] = value.iloc[start:(start + len(data))].values
                    self.data[c] = data
                    start += len(data)
            else:
                raise IndexError('Length of input ({}) does not match data length ({}).'.format(len(value), len(self)))
        if isinstance(value, list):
            if len(value) == len(self):
                start = 0
                for c in range(1, 24):
                    data = self.data[c]
                    data.loc[:, key] = value[start:(start + len(data))]
                    start += len(data)
                    self.data[c] = data
            else:
                raise IndexError('Length of input ({}) does not match data length ({}).'.format(len(value), len(self)))

    def close(self):
        if not self.low_ram:
            print('No open connections')
            return None
        else:
            print("Warning: pickled pysumstats will not work anymore.")
            self.data.close()

    def sort_values(self, by, inplace=True, **kwargs):
        if not inplace:
            raise NotImplementedError()
        if by == 'chr':
            pass
        else:
            for c in self.data.keys():
                data = self.data[c]
                data.sort_values(by, **kwargs)
                self.data[c] = data

    def groupby(self, *args, **kwargs):
        return pd.concat([self.data[c] for c in self.data.keys()]).groupby(*args, **kwargs)

    def save(self, path, per_chromosome=False, **kwargs):
        if ('index' in kwargs.keys()) or ('header' in kwargs.keys()):
            raise KeyError('\'index\' and \'header\' arguments not supported.')
        if not per_chromosome:
            if path.endswith('.pickle'):
                if self.low_ram:
                    warnings.warn(
                        "Saving pysumstats as pickled objects with low_ram will not store the data in the pickled object.")
                with open(path, 'wb') as f:
                    pickle.dump(self, f)
            else:
                with open(path, 'w', newline='', encoding='utf-8') as f:
                    for c, data in self.data.items():
                        if c == 1:
                            data.to_csv(f, index=False, **kwargs)
                        else:
                            data.to_csv(f, header=False, index=False, **kwargs)
        else:
            if (not path.endswith('.csv')) and (not path.endswith('.csv.gz')) and (not path.endswith('.txt')) and (
            not path.endswith('.txt.gz')):
                raise NotImplementedError('Saving files per chromosome only works with .csv(.gz) or .txt(.gz) files')
            for c, data in self.data.items():
                if '{}' not in path:
                    data.to_csv('chr{}_'.format(c) + path, index=False, **kwargs)
                else:
                    data.to_csv(path.format(c), index=False, **kwargs)

    def reset_index(self):
        mn = 0
        mx = 0
        for c in range(1, 24):
            data = self.data[c]
            mx += len(data)
            data.index = list(range(mn, mx))
            mn += len(data)
            self.data[c] = data

    def head(self, n=10, n_chromosomes=1, **kwargs):
        for n_chr in range(1, n_chromosomes + 1):
            print(self.data[n_chr].head(n, **kwargs))

    def tail(self, n=10, n_chromosomes=1, **kwargs):
        for n_chr in range(23, 23 - n_chromosomes):
            print(self.data[n_chr].tail(n, **kwargs))

    def _get_summary(self, sum_cols, per_chromosome):
        summary = {}
        for c in sum_cols:
            summary[c] = self.data[1][c].describe()
        for data in [v for k, v in self.data.items() if k in list(range(2, 24))]:
            for c in summary.keys():
                summary[c] = pd.concat([summary[c], data[c].describe()], axis=1)
        for c in summary.keys():
            summary[c].columns = ['chr{}'.format(x) for x in range(1, 24)]
        if per_chromosome:
            return summary
        else:
            for c, summ in summary.items():
                if 'mean' in summ.index:
                    summ.loc['count', 'full'] = summ.loc['count', :].sum(axis=0)
                    summ.loc['min', 'full'] = summ.loc['min', :].min(axis=0)
                    summ.loc['max', 'full'] = summ.loc['min', :].max(axis=0)
                    n = summ.loc['count', 'full']
                    t = summ.transpose()
                    t['m_w'] = (t['mean'] * t['count']) / n
                    summ.loc['mean', 'full'] = t['m_w'].sum()
                    mn = summ.loc['mean', 'full']
                    t['sd_w'] = (t['std'] ** 2) * (t['count'] - 1) + t['count'] * ((mn - t['mean']) ** 2)
                    summ.loc['std', 'full'] = np.sqrt(t['sd_w'].sum() / (n - 1))
                    summary[c] = summ.loc[['count', 'mean', 'std', 'min', 'max'], 'full']
                else:
                    del summary[c]
            cols = list(summary.keys())
            summary = pd.concat(list(summary.values()), axis=1)
            summary.columns = cols
            return summary
