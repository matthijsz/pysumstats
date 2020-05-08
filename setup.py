from setuptools import setup

setup(name='pysumstats',
      version='0.1',
      description='Package for working with GWAS summary statistics',
      long_description='Python package for reading, combining, meta-analyzing, and saving GWAS summary statistics data.',
      keywords='gwas summary statistics genetics',
      url='https://github.com/matthijsz/pysumstats',
      author='Matthijs D. van der Zee',
      author_email='m.d.vander.zee@vu.nl',
      license='MIT',
      packages=['pysumstats'],
      install_requires=[
          'pandas',
          'tables',
          'numpy',
          'scipy'
      ],
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      include_package_data=True,
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)