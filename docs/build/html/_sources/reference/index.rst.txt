.. _reference_guide:

###############
Reference Guide
###############

Paramagnetic module
-------------------

This module handles the paramagnetic centre by defining the magnetic susceptibility 
tensor and methods for PCS, RDC and PRE calculations.

.. toctree::
    :maxdepth: 1
    
    metal

Protein module
--------------

This module handles the protein structure coordinates and includes methods for
loading a PDB file and calculating atomic properites such as CSA or gyromagnetic ratio

.. toctree::
    :maxdepth: 1
    
    protein

Data I/O module
---------------

This module handles the reading and writing of experimental data.

.. toctree::
    :maxdepth: 1
    
    dataparse

Fitting module
--------------

This module handles the fitting of paramagnetic objects to experimental data.

.. autosummary::
    :toctree: generated
    :template: custom_module_template.rst

    paramagpy.fit
