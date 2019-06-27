#!/usr/bin/env python
import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
	name='paramagpy',
	version='0.4',
	description='Calculate paramagnetic effects in NMR spectra of proteins',
	long_description=long_description,
	long_description_content_type="text/markdown",
	url='https://github.com/henryorton/paramagpy',
	author='Henry Orton',
	author_email='henry.orton@anu.edu.au',
	license='GNU GPLv3+',
	packages=setuptools.find_packages(),
	install_requires=['numpy', 'scipy','matplotlib','biopython'],
	zip_safe=False,
	package_data={'paramagpy': ['paramagpy/icon.gif']},
	include_package_data = True,
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
		"Operating System :: OS Independent"])
