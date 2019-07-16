.. _paramagpy_gui:

Graphic User Interface (GUI)
============================

Paramagpy is equipped with a GUI which is cross-platform and contains most of the functionality of the scripted module. This gives a rapid way for new users to fit and compare PCS, RDC and PRE effects.

YouTube Tutorial
----------------

`Check out the tutorial on YouTube <https://youtu.be/MAoBItSac-g>`_

.. raw:: html

    <iframe width="800" height="450" src="https://www.youtube.com/embed/MAoBItSac-g" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Running the GUI
---------------

To run the GUI, first open the python inperpreter in the terminal

.. code-block:: bash

    user@computer:~$ python3
    Python 3.5.2 (default, Nov 23 2017, 16:37:01) 
    [GCC 5.4.0 20160609] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> 

Then import paramagpy and start the gui with :py:func:`paramagpy.gui.run`.

.. code-block:: bash

    user@computer:~$ python3
    Python 3.5.2 (default, Nov 23 2017, 16:37:01)
    [GCC 5.4.0 20160609] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import paramagpy
    >>> paramagpy.gui.run()


Alternatively you can simply execute the following from the command line

.. code-block:: bash

    user@computer:~$ echo "import paramagpy; paramagpy.gui.run()" | python3

If all this fails, you can contact the author for a prebuilt executable at henry.orton@anu.edu.au

.. image:: paramagpy_gui_annot.png
