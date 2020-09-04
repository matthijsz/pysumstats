"""
.. module:: base
   :synopsis: Base classes used to generate sumstats classes.

.. moduleauthor:: Matthijs van der Zee <m.d.vander.zee@vu.nl>

"""
import pandas as pd
from .exceptions import sumstatswarn
import pickle
import os
import copy
import numpy as np
import time
import gzip
from . import __version__

class _Loc:
    """A helper class to enable sumstats.loc[] similar to pandas dataframes.

    :param master: object this _Loc class is attached to.

    """
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
    """ A helper class to create a connection to a local HDF5 file which can be adressed like a dictionary. Only used locally when low_ram is specified

    :param filename: object this _Loc class is attached to.

    """
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
        """

        :return: range(1, 24); i.e. chromosome numbers

        """
        return range(1, 24)

    def values(self):
        """

        :return: yield data at each chromosome
        """
        for i in self.keys():
            yield self[i]

    def items(self):
        """

        :return: yield (chromsome number, data)
        """
        for i in self.keys():
            yield i, self[i]

    def close(self):
        """Closes connection and removes the file

        """
        os.remove(self.path)


class _BaseSumStats:
    """ Base class of sumstats with functions accessible to both sumstats  and mergedsumstats.

    """
    def __init__(self):
        self.__version__ = __version__
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

    @staticmethod
    def _extract_phenotype(data, phe):
        """Internal function to extract phenotype from MergedSumStats objects.
        :param data: pd.Dataframe
        :param phe: phenotype to extract
        :return: pd.Dataframe with extracted phenotype
        """
        if phe is None:
            return data
        else:
            extracts = ['rsid', 'chr', 'bp', 'ea', 'oa']+[x for x in data.columns if x.endswith('_'+phe)]
            renames = {x: x.replace('_'+phe, '') for x in data.columns if x.endswith('_'+phe)}
            return data[extracts].rename(columns=renames)

    def copy(self):
        """

        :return: a deepcopy of the existing opject
        """
        return copy.deepcopy(self)

    def close(self):
        """Close connection to and HDF5 file if low_ram is specified

        :return: None

        """
        if not self.low_ram:
            print('No open connections')
            return None
        else:
            sumstatswarn("pickled pysumstats saved from this object will not work anymore.")
            self.data.close()

    def sort_values(self, by, inplace=True, **kwargs):
        """Sorts values in the dataframe. Note: Sorting by chromosme (chr) will have no effect as data is already structured by chromosome.

        :param by: label of the column to sort values by
        :type by: str
        :param inplace: Whether to return the sorted object or sort values within existing object. (Currently only inplace sorting is supported)
        :type inplace: bool
        :param kwargs: Other keyword arguments to be passed to pandas sort_values function
        :return: Non

        """
        assert isinstance(by, str), "by should be str"
        assert isinstance(inplace, bool), "inplace should be True or False"
        if not inplace:
            raise NotImplementedError()
        if by == 'chr':
            pass
        elif by not in self.data[1].columns:
            raise AssertionError("{} not found in columns")
        else:
            for c in self.data.keys():
                data = self.data[c]
                data.sort_values(by, **kwargs)
                self.data[c] = data

    def groupby(self, *args, **kwargs):
        """Compatibility function to create pandas grouped object

        :param args: arguments to be passed to pandas groupby function
        :param kwargs: keyword arguments to be passed to pandas groupby function
        :return: a full grouped pandas dataframe object
        """
        return pd.concat([self.data[c] for c in self.data.keys()]).groupby(*args, **kwargs)

    def save(self, path, per_chromosome=False, per_phenotype=False, phenotype=None, **kwargs):
        """Save the data held in this object to local storage.

        :param path: Relative or full path to the target file to store the data or object in. Paths ending in .pickle will save a pickled version of the full object. Note that with low_ram enabled this will **not** store the data. When per_phenotype is specified, add {} to the path where the phenotype name should be, if {} is not in the string, the filename will be prefixed with phenotype name.
        :type path: str
        :param per_chromosome: Whether to save seperate files for each chromosome.
        :type per_chromosome: bool
        :param per_phenotype: Set to True to create a separate file for each phenotype in MergedSumStats objects
        :type per_phenotype
        :param phenotype: Only save a file for a specifici phenotype in MergedSumstats objects
        :type phenotype: str
        :param kwargs: keyword arguments to be passed to pandas to_csv() function.
        :return: None
        """
        assert isinstance(path, str), "path should be str"
        assert isinstance(per_chromosome, bool), "per_chromosome should be True or False"
        if ('index' in kwargs.keys()) or ('header' in kwargs.keys()):
            raise KeyError('\'index\' and \'header\' arguments not supported.')
        if (per_phenotype and (self.phenotype_name is not None)) or ((phenotype is not None) and (self.phenotype_name is not None)):
            sumstatswarn("per_phenotype and phenotype ignored for SumStats object")
            self.save(path, per_chromosome, per_phenotype=False, phenotype=None, **kwargs)
        elif per_phenotype:
            if phenotype is not None:
                sumstatswarn("phenotype ignored when per_phenotype is True")
            if '{}' not in path:
                path = '{}_' + path
            for p in self.pheno_names:
                self.save(path.format(p), per_chromosome, per_phenotype=False, phenotype=p, **kwargs)
        else:
            if not per_chromosome:
                if path.endswith('.pickle'):
                    if self.low_ram:
                        sumstatswarn(
                            "Saving pysumstats as pickled objects with low_ram will not store the data in the pickled object.")
                    with open(path, 'wb') as f:
                        pickle.dump(self, f)
                else:
                    if 'sep' not in kwargs.keys():
                        if (path.endswith('.tsv')) or (path.endswith('.tsv.gz')) or (path.endswith('.txt')) or (path.endswith('.txt.gz')):
                            kwargs['sep'] = '\t'
                        elif (not path.endswith('.csv.gz')) and (not path.endswith('.csv')):
                            raise KeyError('sep should be specified when the requeste output is not .tsv(.gz), .txt(.gz) or .csv(.gz)')
                    if path.endswith('.gz'):
                        with gzip.open(path, 'wb') as f:
                            for c, data in self.data.items():
                                if c == 1:
                                    f.write(self._extract_phenotype(data, phenotype).to_csv(index=False, **kwargs).encode('utf-8'))
                                else:
                                    f.write(self._extract_phenotype(data, phenotype).to_csv(index=False, **kwargs).encode('utf-8'))

                    else:
                        with open(path, 'w', newline='', encoding='utf-8') as f:
                            for c, data in self.data.items():
                                if c == 1:
                                    self._extract_phenotype(data, phenotype).to_csv(f, index=False, **kwargs)
                                else:
                                    self._extract_phenotype(data, phenotype).to_csv(f, header=False, index=False, **kwargs)
            else:
                for c, data in self.data.items():
                    if '{}' not in path:
                        self._extract_phenotype(data, phenotype).to_csv('chr{}_'.format(c) + path, index=False, **kwargs)
                    else:
                        self._extract_phenotype(data, phenotype).to_csv(path.format(c), index=False, **kwargs)

    def reset_index(self):
        """ Reset the index of the data.

        :return: None

        """
        mn = 0
        mx = 0
        for c in range(1, 24):
            data = self.data[c]
            mx += len(data)
            data.index = list(range(mn, mx))
            mn += len(data)
            self.data[c] = data

    def head(self, n=10, n_chromosomes=1, **kwargs):
        """Prints (n_chromosomes) dataframes with the first n rows.

        :param n: number of rows to show
        :type n: int
        :param n_chromosomes: number of chromosomes to show.
        :type n_chromosomes: int
        :param kwargs: keyword arguments to be passed to pandas head function
        :return: None

        """
        assert isinstance(n, int), "n should be int"
        assert isinstance(n_chromosomes, int), "n_chromosomes should be int"
        for n_chr in range(1, n_chromosomes + 1):
            print(self.data[n_chr].head(n, **kwargs))

    def tail(self, n=10, n_chromosomes=1, **kwargs):
        """Prints (n_chromosomes) dataframes with the last n rows.

        :param n: number of rows to show
        :type n: int
        :param n_chromosomes: number of chromosomes to show.
        :type n_chromosomes: int
        :param kwargs: keyword arguments to be passed to pandas tail function
        :return: None

        """
        assert isinstance(n, int), "n should be int"
        assert isinstance(n_chromosomes, int), "n_chromosomes should be int"
        for n_chr in range(23, 23 - n_chromosomes):
            print(self.data[n_chr].tail(n, **kwargs))

    def _get_summary(self, sum_cols, per_chromosome):
        """A function to generate summaries.

        :param sum_cols: columns to generate summaries from
        :type sum_cols: list
        :param per_chromosome: Generate summaries per chromosome
        :type per_chromosome: bool
        :return: pd.DataFrame

        """
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
                    summ.loc['max', 'full'] = summ.loc['max', :].max(axis=0)
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

    def plot_all(self, dest='.', prefix='SumStatsPlots', kwargs={}):
        """Runs all attached plot functions

        :param dest: Folder to save resulting files to. File names will be: {prefix}_{plottype}_{YEAR-MONTH-DAY}.png
        :type dest: str
        :param prefix: prefix to use when saving files.
        :type prefix: str
        :param kwargs: Nested dictionary of other keyword arguments to be passed to each function (keys of top-level dictionary should be function names). Use the 'all' key the top level dictionary to pass keyword argument to every function.
        :type kwargs: dict
        :return: None
        """
        assert isinstance(dest, str), "dest should be str"
        assert isinstance(prefix, str), "prefix should be str"
        assert isinstance(kwargs, dict), "kwargs should be dict"
        if (prefix != '') and (not prefix.endswith('_')):
            prefix += '_'
        if not os.path.isdir(dest):
            raise IOError('Destination is not an existing folder.')
        if 'all' in kwargs.keys():
            for func_key in self.plot_funcs.keys():
                if func_key not in kwargs.keys():
                    kwargs[func_key] = {}
                for k, v in kwargs['all'].items():
                    kwargs[func_key][k] = v
        for funcname, func in self.plot_funcs.items():
            if funcname in kwargs.keys():
                if 'filename' not in kwargs.keys():
                    kwargs['filename'] = '{}/{}{}_{}.png'.format(dest, prefix, funcname, time.strftime('%Y-%m-%d'))
                func(**kwargs[funcname])
        for funcname, func in self.plot_funcs.items():
            if funcname in kwargs.keys():
                if 'filename' not in kwargs[funcname].keys():
                    kwargs[funcname]['filename'] = '{}/{}{}_{}.png'.format(dest, prefix, funcname,
                                                                           time.strftime('%Y-%m-%d'))
                func(**kwargs[funcname])
            else:
                func(filename='{}/{}{}_{}.png'.format(dest, prefix, funcname, time.strftime('%Y-%m-%d')))
