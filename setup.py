#!/usr/bin/env python
from distutils.core import setup

setup(name='paramagpy',
      version='0.1',
      description='Calculate paramagnetic effects in NMR spectra of proteins',
      url='https://github.com/henryorton/paramagpy',
      author='Henry Orton',
      author_email='henry.orton@anu.edu.au',
      license='MIT',
      packages=['paramagpy'],
      requires=['numpy', 'scipy','matplotlib','biopython'],
      zip_safe=False)
