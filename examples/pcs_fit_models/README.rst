.. _pcs_fit_models:

Fit Tensor to PDB with Models
=============================

This example shows how to fit a :math:`{\Delta\chi}`-tensor to experimental PCS data using an NMR structure that contains multiple models. Data for calbindin D9k are used as in the previous example :ref:`pcs_fit`.

There are 3 fitting options available in paramagpy for fitting:

1. Averaged fit: A tensor is fit to each model independently, and then all fitted tensors are averaged together. This is a good choice if models in your PDB represent structural uncertainty.

2. Ensemble averaged fit: A single tensor is fit simultaneously to all models by averaging calculated PCS values during fitting. This is a good choice if models in your PDB represent dynamics as comes from a molecular dynamics simulation.

3. Separate model fit: A tensor is fit to each model independently and the best fitting model is taken. This is a good choice if you are only looking for the best fit model in a PDB containing many models.


Downloads
---------

* Download the data files ``2bcb.pdb`` and ``calbindin_Er_HN_PCS.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit_models.py <../../../examples/pcs_fit_models/pcs_fit_models.py>`


Script + Explanation
--------------------

Firstly, the standard preamble and loading of data.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 1-6

If all models are provided in the ``parsedData`` argument, the default functionality for all fitting methods such as :py:func:`paramagpy.fit.nlr_fit_metal_from_pcs` is to fit using method 1, meaning a tensor is fit to each model and the averaged tensor is returned. This is equivalent to setting the ``ensebleAverage`` argument to ``False``. This is done below. Averaging behaviour can be controlled through the ``idx`` column of ``parsedData``. The ``idx`` array contains common integers for corresponding atoms to be averaged, and defaults to the atom's serial number found in the PDB file.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 12-16

Method 2 can be followed by the same method, except setting the ``ensebleAverage`` argument to ``True``. At each stage of the fitting process, all PCS calculations are then averaged before fitting of a single tensor to all the data simultaneously. The ensemble averaging behaviour can be set through the ``idx`` column of the input data for :py:func:`paramagpy.fit.nlr_fit_metal_from_pcs`.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 18-22

Method 3 can be achieved by constructing a ``for`` loop over the PDB models and fitting a separate tensor to the data from each model. The model which achieves the lowest Q-factor can then be extracted.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 24-31

Finally we plot three sets of data:

    * The averaged fit calculated over all models (green)

    * The ensemble average of the calculated values of the ensemble fit (red)

    * The best fitting single model (blue)

Note that to calculate the ensemble average of the calculated values we use the function :py:func:`paramagpy.fit.ensemble_average`. This can take any number of arguments, and will average values based on common serial numbers of the list of atoms in the first argument.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 34-

*Output:* [:download:`pcs_fit_models.png <../../../examples/pcs_fit_models/pcs_fit_models.png>`]

.. image:: ../../../examples/pcs_fit_models/pcs_fit_models.png