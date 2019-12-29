.. _pcs_fit_models:

Fit Tensor to PDB with Models
=============================

This example shows how to fit a :math:`{\Delta\chi}`-tensor to experimental PCS data using an NMR structure that contains many models. The tensor can be fit to ensemble averaged PCS values, or to individual models. An ensemble averaged PCS is the mean calculated PCS of all models. No structural averages are ever taken.

Data for calbindin D9k are used as in the previous example :ref:`pcs_fit`.


Downloads
---------

* Download the data files ``2bcb.pdb`` and ``calbindin_Er_HN_PCS.npc`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script :download:`pcs_fit_models.py <../../../examples/pcs_fit_models/pcs_fit_models.py>`


Script + Explanation
--------------------

Firstly, the standard preamble and loading of data.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 1-7

The default method of fitting is to fit the tensor independently to each model of the PDB file. To achieve ensemble averaging behaviour, this requires setting the argument ``ensembleAverage`` to true within the fitting function. The default ensemble averaging behaviour is to average atoms with the same serial number in the PDB file. To manipulate ensemble averaging, you can specify the ``idx`` column of the input dataArray for :py:func:`paramagpy.fit.nlr_fit_metal_from_pcs`. The ``idx`` array contains common integers for corresponding atoms to be averaged. After fitting the Q-factor is calculated (with the ensembleAverage argument set to True), and ensemble averaging of the calculated data values is achieved with the function :py:func:`paramagpy.fit.ensemble_average`

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 10-19

Fitting a separate tensor to each model of the PDB is the default behaviour of the fitting functions, the average of all fitted tensors is then returned. The model with the minimum Q-factor is then found by looping over the calculated data and sorting them by the calculated Q-factor.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 22-30

Finally we plot three sets of data:

    * The ensemble average fit calculated for each model (green)

    * The ensemble average of the calculated values of the ensemble fit (red)

    * The best fitting single model (blue)

Note that to calculate the ensemble average of the calculated values we use the function :py:func:`paramagpy.fit.ensemble_average`. This can take any number of arguments, and will average values based on common serial numbers of the list of atoms in the first argument.

.. literalinclude:: ../../../examples/pcs_fit_models/pcs_fit_models.py 
	:lines: 33-

*Output:* [:download:`pcs_fit_models.png <../../../examples/pcs_fit_models/pcs_fit_models.png>`]

.. image:: ../../../examples/pcs_fit_models/pcs_fit_models.png