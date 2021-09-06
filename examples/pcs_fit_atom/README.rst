.. _pcs_fit_atom:

Fit Atomic Coordinates to PCS data
==================================

This example shows how to calculate the region in space which is likely for atomic coordinates from PCS measurements. Multiple :math:`{\Delta\chi}`-tensors are used from different tagging sites in the protein IMP1 to localise a tryptophan sidechain in a loop.


Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Er_HN_PCS.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit.py <../../../examples/pcs_fit/pcs_fit.py>`


Script + Explanation
--------------------

Firstly, the necessary modules are imported from paramagpy. 

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 1

The protein is then loaded from a PDB file using :py:func:`paramagpy.protein.load_pdb` into the variable ``prot``. This returns a ``CustomStructure`` object which is closely based on the ``Structure`` object from `BioPython <https://biopython.org/>`_ and contains the atomic coordinates. The object, and how to access atomic coordinates is discussed at this `link <https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ>`_.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 3-4

The PCS data is then loaded from a ``.npc`` file using the function :py:func:`paramagpy.dataparse.read_pcs` into the variable ``rawData``. This is a dictionary of ``(PCS, Error)`` tuples which may be accessed by ``rawData[(seq, atom)]`` where ``seq`` is an integer specifying the sequence and ``atom`` is the atom name e.g ``(3,'HA')``. Note that these should match the corresponding sequence and atom in the PDB file.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 6-7

To associate the experimental PCS value with atoms of the PDB structure, the method :py:func:`paramagpy.protein.CustomStructure.parse` is called on ``rawData``. The returned array ``parsedData`` has a row for each atom with columns ``[mdl,atm,exp,cal,err,idx]``, where ``mdl`` is the model number from the PDB file, ``atm`` is an atom object from the BioPython PDB structure, ``exp`` and ``cal`` are the experimental and calculated values, ``err`` is the experimental uncertainty and ``idx`` is the atom index, used to define ensemble averaging behaviour.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 9-10

An initial :math:`{\Delta\chi}`-tensor is defined by initialising a :py:class:`paramagpy.metal.Metal` object. The initial position is known to be near the binding site, which is set to the CA atom of residue 56. Note that the ``position`` attribute is always in Angstrom units.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 12-16

A quick gridsearch is conducted in a sphere of 10 Angstrom with 10 points per radius using the function :py:func:`paramagpy.fit.svd_gridsearch_fit_metal_from_pcs`. This requires two lists containing the starting metals ``mStart`` and parsed experimental dataArray ``parsedData``. This function returns lists containing a new fitted metal object, the calculated PCS values from the fitted model.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 18-20

This is then refined using a non-linear regression gradient descent with the function :py:func:`paramagpy.fit.nlr_fit_metal_from_pcs`.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 22-23

The Q-factor is then calculated using the function :py:func`paramagpy.fit.qfactor`.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 25-26

The fitted tensor parameters are saved by calling the method :py:func:`paramagpy.metal.Metal.save`. Alterntaively they may be displayed using ``print(mFit.info())``

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 28-29

*Output:* [:download:`calbindin_Er_HN_PCS_tensor.txt <../../../examples/pcs_fit/calbindin_Er_HN_PCS_tensor.txt>`]

.. literalinclude:: ../../../examples/pcs_fit/calbindin_Er_HN_PCS_tensor.txt
    :language: none

These experimental/calculated PCS values are then plotted in a correlation plot to assess the fit. This is achieved using standard functions of the plotting module `matplotlib <https://matplotlib.org/>`_.

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py 
	:lines: 31-

*Output:* [:download:`pcs_fit.png <../../../examples/pcs_fit/pcs_fit.png>`]

.. image:: ../../../examples/pcs_fit/pcs_fit.png