Opening and saving sumstats files
=================================

Opening sumstats files
----------------------

From csv/txt/tsv
^^^^^^^^^^^^^^^^
pysumstats will automatically recognize csv/txt/tsv extensions and read files appropriately. It assumes a .txt or .tsv files is separated by tabs, and .csv files are separated by commas.
If your file is separated by a character other then the default, please specify `sep` when opening SumStats.
For example, a European csv file can be opened with :class:`pysumstats.SumStats` `('myfile.csv', sep=';')`


From gzipped files
^^^^^^^^^^^^^^^^^^
Since pysumstats uses :func:`pandas.read_csv` to open files the `.gz` extension is automatically recognized and will be opened appropriately.
If you are using a different type of compression see the documentation of :func:`pandas.read_csv` function to check wether your comprression is supported.

From pickled files
^^^^^^^^^^^^^^^^^^
When you have saved your :class:`pysumstats.SumStats` object as a `.pickle` file you can open it as you would any other file.
Loading :class:`pysumstats.MergedSumStats` objects from a `.pickle` file will not work by running `pysumstats.SumStats('mymerged.pickle')`.
To open :class:`pysumstats.MergedSumStats` objects stored as `.pickle` please run the following:

.. code-block:: python

    import pickle
    with open('mymerged.pickle', 'rb') as f:
        my_merged_obj = pickle.load(f)

.. warning::
    When a :class:`pysumstats.SumStats` with `low_ram` enabled is storted as a `.pickle` file the data is **not** included in the file, only the reference to the temporary file!

    When `low_ram` is enabled I highly recommend to save your file as a .csv or something similar.

From other files
^^^^^^^^^^^^^^^^
Any file extension other then the following: `['.txt.', '.txt.gz', '.tsv', '.tsv.gz', '.csv', '.csv.gz', '.pickle']` will be passed straight to :func:`pandas.read_csv` hence you will have to supply the appropriate keyworded arguments to enable reading.

Saving sumstats files
----------------------

Compatibility method
^^^^^^^^^^^^^^^^^^^^^
If you want to use the saved data in any other package or program you should save it in plaintext format (or gzipped plain text).
Both :class:`pysumstats.SumStats` and :class:`pysumstats.MergedSumStats` objects support saving to `.csv, .txt(.gz), .tsv(.gz)` using the same default separators as for reading (commas for .csv, tabs for .txt and .tsv).
Supplying the `sep` argument will override these defaults. Attempting to save with any other file extension will raise a `KeyError`
.. note::

   The `index` and `header` arguments are used internally and cannot be set.


Fastest method
^^^^^^^^^^^^^^^^^^^^^
If you want to continue working on your file with the `pysumstats` package later it is much faster to save your object as a pickled file.
Both :class:`pysumstats.SumStats` and :class:`pysumstats.MergedSumStats` objects support saving to a `.pickle` file from their respective `.save()` methods.
Alternatively files can be pickled in the standard way:

.. code-block:: python

    import pickle
    with open('mypickledsumstats.pickle', 'wb') as f:
        pickle.load(f, my_sumstats_obj)

.. warning::
    When a :class:`pysumstats.SumStats` or :class:`pysumstats.MergedSumStats` with `low_ram` enabled is storted as a `.pickle` file the data is **not** included in the file, only the reference to the temporary file!

    When `low_ram` is enabled I highly recommend to save your file as a .csv or something similar.


Per chromosome
^^^^^^^^^^^^^^
Alternatively, data can be stored in 1 file per chromosome by passing `per_chromsome=True` to :meth:`pysumstats.SumStats.save` or :meth:`pysumstats.SumStats.save`.
By default files will be stored as: `chr1_[path]`, `chr2_[path]`, etc. Alternatively you can add `{}` to the path where you want the chromsomenumbers to be in the files.
For example: `my_sumstats_obj.save('my_data_chr{}.csv', per_chromosome=True)`.