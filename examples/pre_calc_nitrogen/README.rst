.. _pre_calc_nitrogen:

Calculate 15N PREs with cross-correlation effects
=================================================

This example shows how to conduct a weighted fit of a :math:`{\Delta\chi}`-tensor to experimental PCS data with experimental errors.


Downloads
---------

* Download the data files ``4icbH_mut.pdb``, ``calbindin_Tb_N_R1_600.pre`` and ``calbindin_Tb_HN_PCS_tensor.txt`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pre_calc_nitrogen.py <../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.py>`


Script + Explanation
--------------------

First the relevant modules are loaded, the protein and data are read and the data is parsed by the protein.

.. literalinclude:: ../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.py 
	:lines: 1-10

The Tb tensor fitted from PCS data is loaded and the relevant parameters, in this case the magnetic field strength, temperature and rotational correlation time are set.

.. literalinclude:: ../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.py 
	:lines: 13-16

A loop is conducted over the nitrogen atoms that are present in the experimental data. The PRE is calculated using the function :py:func:`paramagpy.metal.atom_pre`. Calculations without CSA are appended to the list ``cal`` and calculations including CSA cross-correlation with the Curie-spin relaxation are appended to the list ``cal_csa``.

.. literalinclude:: ../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.py 
	:lines: 19-25

Finally the data are plotted. Clearly CSA cross-correlation is a big effect for backbone nitrogen atoms and should always be taken into account for Curie-spin calculations. Also note the existence and correct prediction of negative PREs!

.. literalinclude:: ../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.py 
	:lines: 28-45

*Output:* [:download:`pre_calc_nitrogen.png <../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.png>`]

.. image:: ../../../examples/pre_calc_nitrogen/pre_calc_nitrogen.png