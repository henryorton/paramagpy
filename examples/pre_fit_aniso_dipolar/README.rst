.. _pre_fit_aniso_dipolar:

Fit spectral power density tensor
=================================

This example shows how to fit the spectral power density tensor to anisotropic PREs. The data and theory are derived from https://doi.org/10.1039/C8CP01332B.


Downloads
---------

* Download the data files ``parashift_Tb.pdb`` and ``parashift_Tb_R1_exp.pre`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pre_fit_aniso_dipolar.py <../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.py>`


Script + Explanation
--------------------

Load the relevant modules, read the PDB coordinates and experimental PRE values. Parse the values.

.. literalinclude:: ../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.py 
	:lines: 1-7

The spectral power density tensor is written here explicitly and set to the attribute ``g_tensor``. The values here are sourced from the original paper, and arise from the robust linear fit to the experimental data. We will use this tensor for comparison to the fit achieved by paramagpy.

.. literalinclude:: ../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.py 
	:lines: 9-15

An starting tensor with no parameters is also initialised and will be used for fitting to the exerimental data with paramagpy.

.. literalinclude:: ../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.py 
	:lines: 17-18

The fit is conducted by setting the ``usegsbm`` flag to ``True``. This uses anisotropic SBM theory to fit the spectral power density tensor in place of the isotropic SBM theory. The relevant fitting parameters must be specified as ``'t1e', 'gax', 'grh', 'a','b','g'`` which represent the electronic relaxation time, the axial and rhombic componenets of the power spectral density tensor and the 3 Euler angles alpha, beta and gamma respectively. Note that the fitted ``t1e`` parameter is only an estimate of the electronic relaxation time.

.. literalinclude:: ../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.py 
	:lines: 20

Finally the results of the fit are plotted alongside the isotropic theory and the literature fit. Note that the difference in the fit from paramagpy is small, and probably arises because the original paper uses a `Robust` linear fit, which may include weighting with experimental uncertainties. However paramagpy weights values evely here because the experimental uncertainties are unknown.

.. literalinclude:: ../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.py 
	:lines: 24-


*Output:* [:download:`pre_fit_aniso_dipolar.png <../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.png>`]

.. image:: ../../../examples/pre_fit_aniso_dipolar/pre_fit_aniso_dipolar.png