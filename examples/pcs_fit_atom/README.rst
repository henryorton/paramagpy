.. _pcs_fit_atom:

Fit Atomic Coordinates to PCS data
==================================

This example shows how to calculate the region in space which is likely for atomic coordinates from PCS measurements. Multiple :math:`{\Delta\chi}`-tensors are used from different tagging sites in the protein IMP1 to localise a tryptophan sidechain in a loop. 

The script fits the :math:`{\Delta\chi}`-tensors from backbone PCS data and then samples 20 perturbed tensors using a bootstrap fitting. The sampled tensors improve stability of the final calculation. The script then calculates the RMSD between experiment and calculated PCS values for nuclei in the sidechain of tryptophan at residue 28 on a grid of points. The grid is then viewed in PyMol.


Downloads
---------

* Download the data files ``5ev6AH.pdb`` and all IMP1 datasets from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit_atom.py <../../../examples/pcs_fit_atom/pcs_fit_atom.py>`


Script + Explanation
--------------------

After importing modules from paramagpy, the tagging sites and lanthanoid ions are specified as tuple variables. The PDB file is loaded. The ``atoms`` list contains the atoms of interest (in this case the NH and adjacent CH protons of tryptophan 28). Finally the number of bootstrap iterations are defined.

.. literalinclude:: ../../../examples/pcs_fit_atom/pcs_fit_atom.py 
	:lines: 1-10

Two dictionaries are specified to define the final colours and RMSD contour levels to be plotted in PyMol.

.. literalinclude:: ../../../examples/pcs_fit_atom/pcs_fit_atom.py 
	:lines: 12-26

A PyMol script object :py:class:`paramagpy.protein.PyMolScript` is created and the PDB is added to it. This object makes it easy to add density maps, PDBs and spheres to PyMol from Paramagpy.

.. literalinclude:: ../../../examples/pcs_fit_atom/pcs_fit_atom.py 
	:lines: 28-29

Next is a rather involved loop that iterates of the tagging sites, fits the :math:`{\Delta\chi}`-tensor using a simultaneous fit between Tm and Tb data and finally samples the tensor fits using bootstrap. The fitted tensors are bundled into the variable ``mdata``.

.. literalinclude:: ../../../examples/pcs_fit_atom/pcs_fit_atom.py 
	:lines: 31-57

The fitted :math:`{\Delta\chi}`-tensors are then unzipped (to allow iterating over each ion) and assembled with the tryptophan PCS data in two lists ``mdata`` and ``trpdata``. For each data array contained in ``trpdata`` there must be an associated tensor contained in ``mdata``, so that is why they are contsructed side by side.

.. literalinclude:: ../../../examples/pcs_fit_atom/pcs_fit_atom.py 
	:lines: 60-75

The function :py:func:`paramagpy.fit.gridsearch_fit_atom_from_pcs` is called which calculates the PCS RMSD on a grid as defined by the function arguments ``mapSize`` and ``mapDensity``. This function returns a dictionary which contains keys for the atoms of the PDB files and values of :py:class:`paramagpy.fit.DensityMap` which define the grid of PCS RMSD values.

What remains of the script is to add the PCS RMSD grid to the PyMol script and save it so that it plots with the specified colours and contour levels. What results is a volume which contains all points that have an RMSD less than the specified ``isoVals`` value. Finally some standard PyMol commands are added to display the protein structure as desired.

.. literalinclude:: ../../../examples/pcs_fit_atom/pcs_fit_atom.py 
	:lines: 79-

This generates the following PyMol script which allows viewing of the PCS RMSD region. :download:`pcs_fit_atom.pml <../../../examples/pcs_fit_atom/pcs_fit_atom.pml>`. After opening the script in PyMol the following image is generated.

[:download:`pcs_fit_atom.png <../../../examples/pcs_fit_atom/pcs_fit_atom.png>`]

.. image:: ../../../examples/pcs_fit_atom/pcs_fit_atom.png
