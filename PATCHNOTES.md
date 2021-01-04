# Patch notes

##### 04-01-2021 (v0.5)
 - Please update to this version! Merging and meta analyses did not work properly in previous versions.
 - For now using `MergedSumStats.meta_analyze()` will instead run `.gwama()` with an identity matrix (functionally identical to `.meta_analyze(method='ivw')`), as this appears to work reliably.

##### 30-11-2020 (v0.4.6)
 - Fix: Merge using a list of SumStats now propagates arguments properly
 - Fix: Meta analyzing now properly sums N across all summary statistics
 - Fix: Fixed an issue that would prevent proper outer merge

##### 30-11-2020 (v0.4.5)
 - Fix: Merge using a list of SumStats now propagates arguments properly
 - Fix: Meta analyzing now properly sums N across all summary statistics

##### 16-09-2020 (v0.4.4)
 - Bugfix: Chromosome X now properly converted to 23

##### 04-09-2020 (v0.4.3)
 - Bugfix: .describe() now returns min and max instead of min and min
 - Bugfix: .qc() now works with any column
 
##### 28-07-2020 (v0.4.2)
 - Added an NotImplementedError when attempting to merge with `low_memory` and a `method` other than inner.
 - Added per_phenotype argument to `.save()` to save separate files for each phenotype in `MergedSumStats` objects.
 - Added phenotype argument to `.save()` to save file for a specific phenotype in `MergedSumstats` objects.

##### 02-06-2020 (v0.4.1)
 - Fixed assertion issues that caused merging files and plotting merged files to fail.

##### 15-05-2020 (v0.4)
 - Fixed export so that saving to a file with '.gz' extension now actaully gzips the output file.
 - Fixed data import to enable opening a pysumstats file saved as .pickle by initializing SumStats class. This will also check if the data is still stored when low_ram was set to True.
 - Merging sumstats now computes an overall MAF meta-analyzed sumstats can be reimported with pysumstats.
 - Added version number to SumStats objects.
 - Added `kwargs` to SumStats.qc for filtering on custom columns.
 - Added custom SumStats warnings. These warnings print always by default you can change this behavior with Pythons default [warnings](https://docs.python.org/3/library/warnings.html) filter.
 - Reformatted [documentation](https://pysumstats.readthedocs.io/en/latest/), now also includes a guide on opening and saving sumstats.

##### 13-05-2020 (v0.3.1)
 - Fixed an issue where reading data would fail when values in n, bp, chr columns were NA. An attempt is now made to impute these values. If too many are missing a ValueError is thrown.

##### 12-05-2020 (v0.3)
 - Added `fig` and `ax` arguments to `pysumstats.plot.qqplot` and `pysumstats.plot.manhattan` to enable plotting to existing figure and axis.
 - Added `pysumstats.plot.pzplot`, to visually compare Z-values from `B/SE` to Z-values calculated from the P-value.
 - Added `pysumstats.plot.afplot`, to plot allele frequency differences between summary statistics.
 - Added `pysumstats.plot.zzplot`, to plot differences in Z-values between summary statistics.
 - Added `qqplot`, `manhattan`, `pzplot`, `afplot`, `zzplot` functions to MergedSumStats object.
 - Added `pzplot` function to SumStats object.
 - Added `plot_all` functions to SumStats and MergedSumStats objects to automatically generate all possible plots for the object.

##### 11-05-2020 (v0.2.3)
 - Fixed import errors and added `manhattan` and `qq` function to `SumStats` class
 - Added `return` statement to MergedSumStats.merge() when `inplace=False` and merging with other MergedSumstats.
 - Added docstrings to base, mergedsumstats, sumstats and utils.
 - Added [docs](https://pysumstats.readthedocs.io/en/latest/)


##### 08-05-2020 (v0.2)

 - Added `plot` subpackage with `qqplot` and `manhattan`,  from  my initial [Python-QQMan module](https://github.com/matthijsz/qqman).

##### 08-05-2020 (v0.1)

 - Adapted to be a package rather then a module.
 - Added `low_ram` argument to SumStats to read/write data to disk rather than RAM, in case of memory issues.  
