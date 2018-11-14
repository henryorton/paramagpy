.. _pcs_fit_multiple:

Fit multiple PCS datasets to common position
============================================

This example shows how to fit multiple :math:`{\Delta\chi}`-tensors to their respective datasets with a common position, but varied magnitude and orientation. This may arise if several lanthanides were investigated at the same binding site, and the data may be used simultaneously to fit a common position. Data from several PCS datasets for calbindin D9k were used here, and is a generalisation of the previous example: :ref:`pcs_fit`.


Downloads
---------

* Download the data files ``4icbH_mut.pdb``, ``calbindin_Tb_HN_PCS.npc``, ``calbindin_Er_HN_PCS.npc`` and ``calbindin_Yb_HN_PCS_tensor.txt`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script `pcs_fit_multiple.py <https://github.com/henryorton/paramagpy/tree/master/examples/pcs_fit_multiple/pcs_fit_multiple.py>`_


Explanation
-----------

The protein and PCS datasets are loaded and parsed. These are placed into a list ``parsedData``, for which each element is a PCS dataset of a given lanthanide.

The two fitting functions: 

* :py:func:`paramagpy.fit.svd_gridsearch_fit_metal_from_pcs`  

* :py:func:`paramagpy.fit.nlr_fit_metal_from_pcs` 

can accept a list of metal objects and a list of datasets with arbitrary size. If this list contains more than one element, fitting will be performed to a common position. The starting position is taken only from the first metal of the list.

After fitting, a list of fitted metals is returned. The fitted tensor are then written to files and a correlation plot is made.


Script
------

[:download:`pcs_fit_multiple.py <../../../examples/pcs_fit_multiple/pcs_fit_multiple.py>`]

.. literalinclude:: ../../../examples/pcs_fit_multiple/pcs_fit_multiple.py


Outputs
-------

*Tb fitted tensor*

[:download:`tensor_Tb.txt <../../../examples/pcs_fit_multiple/tensor_Tb.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_multiple/tensor_Tb.txt

*Er fitted tensor*

[:download:`tensor_Er.txt <../../../examples/pcs_fit_multiple/tensor_Er.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_multiple/tensor_Er.txt

*Yb fitted tensor*

[:download:`tensor_Yb.txt <../../../examples/pcs_fit_multiple/tensor_Yb.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_multiple/tensor_Yb.txt


*Correlation Plot*

[:download:`pcs_fit_multiple.png <../../../examples/pcs_fit_multiple/pcs_fit_multiple.png>`]

.. image:: ../../../examples/pcs_fit_multiple/pcs_fit_multiple.png