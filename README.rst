paramagpy 
=========

paramagpy is a python module for calculating paramagnetic effects in NMR spectra of proteins. This currently includes fitting of paramagnetic susceptibility tensors to experimental data associated with pseudocontact shifts (PCS) residual dipolar couplings (RDC) and paramagnetic relaxation enhancements (PRE). A GUI allows easy viewing of data and seamless transition between PCS/RDC/PRE calculations.

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

Documentation
-------------

* https://henryorton.github.io/paramagpy/


Citing paramagpy
----------------

Hopefully paramagpy will be published soon the journal of biomolecular NMR