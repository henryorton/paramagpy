.. _fit_pcs_single:

Fit Tensor to PCS Data
======================

This example shows how to fit a :math:`{\Delta\Chi}` tensor to experimental PCS data for the protein calbindin D9k. 

First the protein is loaded from a PDB file using :py:func:`paramagpy.protein.load_pdb`


Instructions
------------

* Download the data files `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files>`_.


``fit_pcs.py`` [:download:`source code <../../../examples/pcs_single_fit/fit_pcs.py>`]

.. literalinclude:: ../../../examples/pcs_single_fit/fit_pcs.py

Output:

[:download:`pcs_plot.pdf <../../../examples/pcs_single_fit/pcs_plot.pdf>`]

.. image:: ../../../examples/pcs_single_fit/pcs_plot.png