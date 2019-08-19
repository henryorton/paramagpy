.. _ccr_calculate:

Calculate Cross-correlated Relaxation
=====================================

This example shows how to calculate dipole-dipole/Curie-spin cross-correlated relaxation as measured for data in the literature by `Pintacuda et. al. <https://doi.org/10.1023/A:1024926126239>`_


Downloads
---------

* Download the data files ``1bzrH.pdb``, ``myoglobin_cn.ccr`` and ``myoglobin_f.ccr`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`ccr_calculate.py <../../../examples/ccr_calculate/ccr_calculate.py>`


Script + Explanation
--------------------

First the relevant modules are loaded, and the iron atom (paramagnetic centre) is identified as the variable ``ironAtom``.

.. literalinclude:: ../../../examples/ccr_calculate/ccr_calculate.py 
	:lines: 1-6

Two paramagnetic centres are defined for the high and low spin iron atom. The positions are set to that of the iron centre along with other relevant parameters. The measured isotropic :math:`{\chi}`-tensor magnitudes are also set.

.. literalinclude:: ../../../examples/ccr_calculate/ccr_calculate.py 
	:lines: 9-15

The experimental data are loaded and parsed by the protein.

.. literalinclude:: ../../../examples/ccr_calculate/ccr_calculate.py 
	:lines: 18-19

A loop is conducted over the atoms contained in the experimental data and the CCR rate is calculated using the function :py:meth:`paramagpy.metal.Metal.atom_ccr`. These are appended to lists ``compare_cn`` and ``compare_f``.

Note that the two H and N atoms are provided. The first atom is the nuclear spin undergoing active relaxation. The second atom is the coupling partner. Thus by swapping the H and N atoms to give ``atom_ccr(N, H)``, the differential line broadening can be calculated in the indirect dimension.

.. literalinclude:: ../../../examples/ccr_calculate/ccr_calculate.py 
	:lines: 21-30

Finally a correlation plot is made.

.. literalinclude:: ../../../examples/ccr_calculate/ccr_calculate.py 
	:lines: 32-

*Output:* [:download:`ccr_calculate.png <../../../examples/ccr_calculate/ccr_calculate.png>`]

.. image:: ../../../examples/ccr_calculate/ccr_calculate.png