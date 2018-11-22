.. _pcs_fit_constrained:

Constrained Fitting
===================

This example shows how to fit a :math:`{\Delta\chi}`-tensor with constraints applied. The two cases here constrain position to fit a tensor to a known metal ion position form an X-ray structure, and fit an axially symmetric tensor with only 6 of the usual 8 parameters.


Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Er_HN_PCS.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit_constrained.py <../../../examples/pcs_fit_constrained/pcs_fit_constrained.py>`


Script + Explanation
--------------------

The necessary modules are imported and data is loaded 

.. literalinclude:: ../../../examples/pcs_fit_constrained/pcs_fit_constrained.py 
	:lines: 1-7

The calcium ion from the X-ray structure is contained in a heteroatom of the PDB file. We set the starting position of the tensor to this position.

.. literalinclude:: ../../../examples/pcs_fit_constrained/pcs_fit_constrained.py 
	:lines: 9-10

To fit the the anisotropy and orientation without position, the linear PCS equation can be solved analytically by the SVD gridsearch method but using only one point with a radius of zero. This tensor is then saved.

.. literalinclude:: ../../../examples/pcs_fit_constrained/pcs_fit_constrained.py 
	:lines: 12-16

*Output:* [:download:`pcs_fit_constrained.png <../../../examples/pcs_fit_constrained/calbindin_Er_HN_PCS_tensor_position_constrained.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_constrained/calbindin_Er_HN_PCS_tensor_position_constrained.txt

To fit an axially symmetric tensor, we can used the Non-linear regression method and specify exactly which parameters we want to fit. This will be the axiality ``ax``, two Euler angles ``b`` and ``g`` and the position coordinates. Note that in the output, the rhombic ``rh`` and alpha ``a`` parameters are redundant.

.. literalinclude:: ../../../examples/pcs_fit_constrained/pcs_fit_constrained.py 
	:lines: 19-23

*Output:* [:download:`pcs_fit_constrained.png <../../../examples/pcs_fit_constrained/calbindin_Er_HN_PCS_tensor_axially_symmetric.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_constrained/calbindin_Er_HN_PCS_tensor_axially_symmetric.txt

Finally we plot the data.

.. literalinclude:: ../../../examples/pcs_fit_constrained/pcs_fit_constrained.py 
	:lines: 26-

*Output:* [:download:`pcs_fit_constrained.png <../../../examples/pcs_fit_constrained/pcs_fit_constrained.png>`]

.. image:: ../../../examples/pcs_fit_constrained/pcs_fit_constrained.png