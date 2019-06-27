.. _pre_fit_proton:


Fit Tensor to PRE Data
======================

This example demonstrates fitting of the rotational correlation time :math:`{\tau_r}` to 1H PRE data of calbindin D9k. You can fit any parameters of the :math:`{\chi}`-tensor you desire, such as position or magnitude as well.


Downloads
---------

* Download the data files ``4icbH_mut.pdb``, ``calbindin_Er_H_R2_600.npc`` and ``calbindin_Tb_H_R2_800.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pre_fit_proton.py <../../../examples/pre_fit_proton/pre_fit_proton.py>`


Script + Explanation
--------------------

Firstly, the necessary modules are imported from paramagpy. 

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 1

The protein is then loaded from a PDB file.

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 3-4

The PRE data is loaded. Note that the Er data was recorded at 600 MHz and the Tb data was recorded at 800 MHz.

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 7-8

The :math:`{\Delta\chi}`-tensors that were fitted from PCS data are loaded from file and the relevant :math:`{B_0}` magnetic field strengths are set.

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 15-18

Fitting of the rotational correlation time is done with the function :py:func:`paramagpy.fit.nlr_fit_metal_from_pre`. To fit position or :math:`{\chi}`-tensor magnitude, you can change the ``params`` argument.

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 21-24

The fitted tensors are saved to file. Note that the Er dataset gives a reasonable :math:`{\tau_r}` of around 4 ns which is close to the literature value of 4.25 ns. However, the Tb dataset gives an unreasonably large value of 18 ns. This is due to magnetisation attenuation due to 1H-1H RDCs present during the relaxation evolution time as discussed in `literature <https://doi.org/10.1021/jacs.8b03858>`_ giving rise to artificially large measured PREs for lanthanides with highly anisotropic :math:`{\Delta\chi}`-tensors. This is also reflected in the correlation plot below.

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 27-28

*Output:* [:download:`calbindin_Er_H_R2_600_tensor.txt <../../../examples/pre_fit_proton/calbindin_Er_H_R2_600_tensor.txt>`]

.. literalinclude:: ../../../examples/pre_fit_proton/calbindin_Er_H_R2_600_tensor.txt
    :language: none

*Output:* [:download:`calbindin_Tb_H_R2_800_tensor.txt <../../../examples/pre_fit_proton/calbindin_Tb_H_R2_800_tensor.txt>`]

.. literalinclude:: ../../../examples/pre_fit_proton/calbindin_Tb_H_R2_800_tensor.txt
    :language: none

And the results are plotted.

.. literalinclude:: ../../../examples/pre_fit_proton/pre_fit_proton.py 
	:lines: 31-54

.. image:: ../../../examples/pre_fit_proton/pre_fit_proton.png