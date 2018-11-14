.. _pcs_fit:

Fit Tensor to PCS Data
======================

This example shows how to fit a :math:`{\Delta\chi}`-tensor to experimental PCS data for the protein calbindin D9k. These data contain amide 1H and 15N chemical shifts between diamagnetic and paramagnetic states with the lanthanide Er3+ bound.


Downloads
---------

* Download the data files ``4icbH_mut.pdb`` and ``calbindin_Er_HN_PCS.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script `pcs_fit.py <https://github.com/henryorton/paramagpy/tree/master/examples/pcs_fit/pcs_fit.py>`_


Explanation
-----------

First the protein is loaded from a PDB file using :py:func:`paramagpy.protein.load_pdb` into the variable ``prot``. This returns a ``CustomStructure`` object which is closely based on the ``Structure`` object from `BioPython <https://biopython.org/>`_ and contains the atomic coordinates. The object, and how to access atomic coordinates is discussed at this `link <https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ>`_.

The PCS data is then loaded from a ``.npc`` file using the function :py:func:`paramagpy.dataparse.read_pcs` into the variable ``rawData``. This is a dictionary of ``(PCS, Error)`` tuples which may be accessed by ``rawData[(seq, atom)]`` where ``seq`` is an integer specifying the sequence and ``atom`` is the atom name e.g ``(3,'HA')``. Note that these should match the corresponding sequence and atom in the PDB file.

To associate the experimental PCS value with atoms of the PDB structure, the method :py:func:`paramagpy.protein.CustomStructure.parse` is called on ``rawData``. The new list ``parsedData`` contains elements ``[atom, PCS, Error]``, where ``atom`` is now an atom object from the PDB.

An initial :math:`{\Delta\chi}`-tensor is defined by initialising a :py:class:`paramagpy.metal.Metal` object. The initial position is known to be near the binding site, which is set to the CA atom of residue 56. Note that the ``position`` attribute is always in Angstrom units.

A quick gridsearch is conducted in a sphere of 10 Angstrom with 10 points per radius using the function :py:func:`paramagpy.fit.svd_gridsearch_fit_metal_from_pcs`. This requires two lists containing the starting metals ``mStart`` and parsed experimental data ``parsedData`` and return a new metal object which fits best.

This is then refined using a non-linear regression gradient descent with the function :py:func:`paramagpy.fit.nlr_fit_metal_from_pcs`.

The fitted tensor parameters are then output by calling the method :py:func:`paramagpy.metal.Metal.save`.


Script
------

[:download:`pcs_fit.py <../../../examples/pcs_fit/pcs_fit.py>`]

.. literalinclude:: ../../../examples/pcs_fit/pcs_fit.py

Output
------

*The fitted tensor parameters*

[:download:`calbindin_Er_HN_PCS_tensor.txt <../../../examples/pcs_fit/calbindin_Er_HN_PCS_tensor.txt>`]

.. literalinclude:: ../../../examples/pcs_fit/calbindin_Er_HN_PCS_tensor.txt