paramagpy 
=========

paramagpy is a python module for calculating paramagnetic effects in NMR spectra of proteins. This currently includes fitting of paramagnetic susceptibility tensors to experimental data associated with pseudocontact shifts (PCS) residual dipolar couplings (RDC), paramagnetic relaxation enhancements (PRE) and cross-correlated relaxation (CCR). A GUI allows easy viewing of data and seamless transition between PCS/RDC/PRE/CCR calculations.

.. figure:: paramagpy/icon.png
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
* Error analysis of tensor fit quality by Monte-Carlo or Bootstrap methods
* Optimisation of multiple PCS/PRE/CCR datasets to a common position
* Unique tensor representation compatible with Numbat (program)
* Fitting of RDC tensor by SVD algorithm
* PRE calculations by Solomon and Curie spin mechanisms
* Spectral power density tensor fitting for anisotropic dipolar PREs
* CSA cross-correlation correction to PRE calculations
* Dipole-dipole/Curie spin cross-correlated relaxation calculations
* Fitting of tensor parameters to PRE/CCR data
* Macro scripts for integration with CCPNMR and Sparky

Documentation
-------------

* https://henryorton.github.io/paramagpy/


Citing paramagpy
----------------

Paramagpy is published in Magnetic Resonance https://doi.org/10.5194/mr-2019-3
