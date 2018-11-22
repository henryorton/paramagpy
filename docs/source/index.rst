.. paramagpy documentation master file, created by
   sphinx-quickstart on Fri Nov  9 13:20:36 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to paramagpy's documentation
====================================

:Release: |version|
:Date: |today|

Contents
--------

.. toctree::
    :maxdepth: 2

    install
    examples/index
    paramagpy_gui
    reference/index

Introduction
------------

Paramagpy is a python module for calculating paramagnetic effects in NMR spectra of proteins. This currently includes fitting of paramagnetic susceptibility tensors to experimental data associated with pseudocontact shifts (PCS) residual dipolar couplings (RDC) and paramagnetic relaxation enhancements (PRE). A GUI allows easy viewing of data and seamless transition between PCS/RDC/PRE calculations.

.. figure:: ../../paramagpy/icon.gif
    :align: center

    *Please, not the eyes!* - Canberra cyclist

Features
--------

* Support for PDB protein structures with models
* Combined SVD gridsearch and gradient descent algorithms for solving PCS tensors
* Optional fitting of reference offset parameter for PCS datasets
* Support for Residual Anisotropic Chemical Shielding (RACS) and Residual Anisotropic Dipolar Shielding (RADS) corrections to PCS
* Lanthanide parameter templates available
* Plotting of correlation between experiment/calculated values
* Plotting of tensor isosurfaces compatible with PyMol
* Q-factor calculations
* Optimisation of multiple PCS/PRE datasets to a common position
* Unique tensor representation compatible with Numbat (program)
* Fitting of RDC tensor by SVD algorithm
* PRE calculations by Solomon and Curie spin mechanisms
* CSA cross-correlation correction to PRE calculations

