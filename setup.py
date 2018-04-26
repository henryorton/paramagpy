from setuptools import setup

setup(name='pylanth',
      version='0.0.1',
      description='Calculate paramagnetic effects in NMR spectra of proteins',
      url='https://github.com/henryorton/pylanth',
      author='Henry Orton',
      author_email='henry.orton@anu.edu.au',
      license='MIT',
      packages=['numpy','scipy','matplotlib'],
      zip_safe=False)