.. _pcs_calc:

Calculate PCSs from known tensor
================================

This example shows how calculate PCSs from a known :math:`{\Delta\chi}`-tensor that is stored in a file.

Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Tb_HN_PCS_tensor.txt`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_calc.py <../../../examples/pcs_calc/pcs_calc.py>`


Script + Explanation
--------------------

This simple script reads the :math:`{\Delta\chi}`-tensor from a file and calcualtes the PCS for all atoms in a PDB file. The calculated PCS is then written to an ``.npc`` file.

.. literalinclude:: ../../../examples/pcs_calc/pcs_calc.py 
	:lines: 1-
