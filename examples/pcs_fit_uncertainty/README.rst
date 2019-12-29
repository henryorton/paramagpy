.. _pcs_fit_uncertainty:

Propagate Uncertainty to Fitted Tensor Parameters
=================================================

This example shows the various error analysis functions available in paramagpy for estimating the unceratinty in fitted parameters for a paramagnetic center.


Downloads
---------

* Download the data files ``2bcb.pdb`` and ``calbindin_Er_HN_PCS_errors.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit_uncertainty.py <../../../examples/pcs_fit_uncertainty/pcs_fit_uncertainty.py>`


Script + Explanation
--------------------

This start of this script follows the script :ref:`pcs_fit` to fit the tensor.

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/pcs_fit_uncertainty.py
	:lines: 1-24

### Uncertainty from structure models

The PDB file contains models that capture uncertainty in the structure of the protein. This can be propagated to estimate uncertainty in the fitted tensor parameters using the fnction :py:func:`paramagpy.fit.fit_error_model`. This fits a separate tensor to each model and returns all fitted tensors as well as the standard deviation in the fitted parameters.

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/pcs_fit_uncertainty.py
	:lines: 26-30

The standard deviation in the fitted tensor parameters is found in the variable ``mod_std``. This variation in tensor principle axes can be viewed by a Sanson-Flamsteed plot.

*Output:* [:download:`error_tensor_models.txt <../../../examples/pcs_fit_uncertainty/error_tensor_models.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/error_tensor_models.txt

*Output:* [:download:`models.png <../../../examples/pcs_fit_uncertainty/models.png>`]

.. image:: ../../../examples/pcs_fit_uncertainty/models.png


### Uncertainty from experimental uncertainties

Experimental uncertainties can be measured. This may arise due to spectral noise in peak heights for PREs, or spectral noise as uncertainties in chemical shifts for PCSs, as is the case here. The function :py:func:`paramagpy.fit.fit_error_monte_carlo` will repeat the fit for many iterations, each time adding random noise from a uniform distribution scaled by the experimental errors present in the ``err`` column of the dataArray ``parsedData``. 

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/pcs_fit_uncertainty.py
	:lines: 32-36

*Output:* [:download:`error_tensor_monte_carlo.txt <../../../examples/pcs_fit_uncertainty/error_tensor_monte_carlo.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/error_tensor_monte_carlo.txt

*Output:* [:download:`monte_carlo.png <../../../examples/pcs_fit_uncertainty/monte_carlo.png>`]

.. image:: ../../../examples/pcs_fit_uncertainty/monte_carlo.png


### Uncertainty from sample fraction

A final, but generally not recommended method is to source noise from taking a random fraction of the data and conducting the fit for many iterations to then view the deviation in fitted parameters. This method is often called bootstrapping and is desirable if the experimental uncertainties are unknown and the PDB file does not contain models that capture structural unceratinty. The function :py:func:`paramagpy.fit.fit_error_bootstrap` will repeat the fit for many iterations, each time sampling the desired amount of the experimental data randomly. 

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/pcs_fit_uncertainty.py
	:lines: 38-42

*Output:* [:download:`error_tensor_bootstrap.txt <../../../examples/pcs_fit_uncertainty/error_tensor_bootstrap.txt>`]

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/error_tensor_bootstrap.txt

*Output:* [:download:`bootstrap.png <../../../examples/pcs_fit_uncertainty/bootstrap.png>`]

.. image:: ../../../examples/pcs_fit_uncertainty/bootstrap.png

This piece of code is used to generate the Sanson-Flamsteed projection plots

.. literalinclude:: ../../../examples/pcs_fit_uncertainty/pcs_fit_uncertainty.py
	:lines: 44-



















