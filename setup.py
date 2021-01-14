from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'README.md'), 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(name='pysumstats',
      version='0.5.1',
      description='Package for working with GWAS summary statistics',
      long_description=long_description,
      long_description_content_type='text/markdown',
      keywords='gwas summary statistics genetics',
      url='https://github.com/matthijsz/pysumstats',
      author='Matthijs D. van der Zee',
      author_email='m.d.vander.zee@vu.nl',
      license='MIT',
      packages=['pysumstats', 'pysumstats.plot'],
      install_requires=[
          'pandas',
          'tables',
          'numpy',
          'scipy',
          'matplotlib'
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
