.. _rdc_calculate:

Calculate RDC from a known Tensor
=================================

This example shows how to calculate theoretical RDC values from a known :math:`{\Delta\chi}`-tensor which has been fitted from PCS data. Paramagpy allows seamless calculation of one PCS/PRE/RDC/CCR effect from a tensor fitted from another effect.


Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Er_HN_PCS_tensor.txt`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`rdc_calculate.py <../../../examples/rdc_calculate/rdc_calculate.py>`


Script + Explanation
--------------------

First the relevant modules are loaded, the protein is loaded and the metal is loaded from file. The magnetic field strength and temperature are also set.

.. literalinclude:: ../../../examples/rdc_calculate/rdc_calculate.py
	:lines: 1-8

A loop is made over the atoms of the protein. The amide H and N atoms are selected and then the RDC value is calculated. Finally the formated data is appended to list ``forFile``.

.. literalinclude:: ../../../examples/rdc_calculate/rdc_calculate.py
	:lines: 12-23

The formatted data is written to file:

.. literalinclude:: ../../../examples/rdc_calculate/rdc_calculate.py
	:lines: 26-27

*Output:* [:download:`calbindin_Er_RDC_calc.rdc <../../../examples/rdc_calculate/calbindin_Er_RDC_calc.rdc>`]

.. literalinclude:: ../../../examples/rdc_calculate/calbindin_Er_RDC_calc.rdc