.. _pcs_fit_error:

Fit a tensor to PCS data with uncertainties
===========================================

This example shows how to conduct a weighted fit of a :math:`{\Delta\chi}`-tensor to experimental PCS data with experimental errors.


Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Er_HN_PCS_errors.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit_error.py <../../../examples/pcs_fit_error/pcs_fit_error.py>`


Script + Explanation
--------------------

This script follows very closely the script :ref:`pcs_fit`. The only difference being that errors are included in the fourth column of the .npc file and errorbars are included in the plotting routine.

.. literalinclude:: ../../../examples/pcs_fit_error/pcs_fit_error.py 

The fitted tensor:

*Output:* [:download:`calbindin_Er_HN_PCS_tensor_errors.txt <../../../examples/pcs_fit_error/calbindin_Er_HN_PCS_tensor_errors.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_error/calbindin_Er_HN_PCS_tensor_errors.txt

And correlation plot:

*Output:* [:download:`pcs_fit_error.png <../../../examples/pcs_fit_error/pcs_fit_error.png>`]

.. image:: ../../../examples/pcs_fit_error/pcs_fit_error.png