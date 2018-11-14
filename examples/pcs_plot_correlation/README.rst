.. _pcs_plot_correlation:

Plot experiment/calculated PCS and Q-factor
===========================================

This example shows how to plot calculated vs experimental PCS values of a fitted :math:`{\Delta\chi}`-tensor and calculate the Q-factor for data from the example :ref:`pcs_fit`.


Downloads
---------

* Download the data files ``4icbH_mut.pdb``, ``calbindin_Er_HN_PCS.npc`` and ``calbindin_Er_HN_PCS_tensor.txt`` from `here <https://github.com/henryorton/paramagpy/tree/master/examples/data_files/>`_:

* Download the script `pcs_plot_correlation.py <https://github.com/henryorton/paramagpy/tree/master/examples/pcs_plot_correlation/pcs_plot_correlation.py>`_


Explanation
-----------

The protein and data are loaded as described previously in :ref:`pcs_fit`.

The fitted metal is loaded into the variable ``met`` from a file using the function :py:func:`paramagpy.metal.load_tensor`.

A *for loop* is used to make two lists, ``exp`` and ``cal`` which contain the experimental and calculated PCS values respectively.

The Q-factor is then calculated using the function :py:func:`paramagpy.fit.qfactor`.

Finally the results are plotted using standard functions of the plotting module `matplotlib <https://matplotlib.org/>`_


Script
------

[:download:`pcs_plot_correlation.py <../../../examples/pcs_plot_correlation/pcs_plot_correlation.py>`]

.. literalinclude:: ../../../examples/pcs_plot_correlation/pcs_plot_correlation.py


Output
------

*Correlation Plot*

[:download:`pcs_plot_correlation.png <../../../examples/pcs_plot_correlation/pcs_plot_correlation.png>`]

.. image:: ../../../examples/pcs_plot_correlation/pcs_plot_correlation.png