.. _pcs_plot_isosurface:

Plot PCS isosurface (PyMol view)
================================

This example shows how to plot the PCS isosurface of a fitted :math:`{\Delta\chi}`-tensor for data from the example :ref:`pcs_fit`. The isosurface can be viewed in `PyMol <https://pymol.org>`_.


Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Er_HN_PCS_tensor.txt`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script `pcs_plot_isosurface.py <https://github.com/henryorton/paramagpy/tree/master/examples/pcs_plot_isosurface/pcs_plot_isosurface.py>`_


Explanation
-----------

The protein and tensor are loaded as described previously in.

The isosurface files are generated using the function :py:func:`paramagpy.metal.Metal.isomap`. The contour level can be chosen by setting the ``isoval`` argument. A larger ``density`` value will result in a smoother surface. This function writes two files ``isomap.pml`` and ``isomap.pml.ccp4`` which are the PyMol script and PCS grid files respectively.

The isosurface can be displayed by executing ``pymol isomap.pml`` from a terminal, or by selecting ``File>Run`` and navigating to the script ``isomap.pml``.


Script
------

[:download:`pcs_plot_isosurface.py <../../../examples/pcs_plot_isosurface/pcs_plot_isosurface.py>`]

.. literalinclude:: ../../../examples/pcs_plot_isosurface/pcs_plot_isosurface.py


Output
------

*PyMol view of isosurface*

[:download:`pcs_plot_isosurface.png <../../../examples/pcs_plot_isosurface/pcs_plot_isosurface.png>`]

.. image:: ../../../examples/pcs_plot_isosurface/pcs_plot_isosurface.png
