.. solsticepy documentation master file, created by
   sphinx-quickstart on Fri Apr 17 13:39:21 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

solsticepy
==========

`solsticepy` provides a set of Python functions that make the task of setting
up and running a Solstice simulation of a CSP system a little easier. At this
stage it has been primarily used for simulations of central-tower CSP systems, 
even though Solstice itself is capable of simulating a wider range of system
types.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Create the YAML input file for Solstice
---------------------------------------

.. autofunction:: solsticepy.gen_yaml

Calculation sun position
------------------------

.. autoclass:: solsticepy.SunPosition
   :members:
   :undoc-members:

Set up and run Solstice simulations
-----------------------------------

.. autoclass:: solsticepy.Master
   :members:
   :undoc-members:

Process the results from Solstice
---------------------------------

.. autofunction:: solsticepy.process_raw_results

Generate 3D views of your system
--------------------------------

.. autofunction:: solsticepy.gen_vtk

Generate new heliostat field layouts
------------------------------------

.. autofunction:: solsticepy.radial_stagger

Preliminary calculation of heliostat field performance
------------------------------------------------------

.. autoclass:: solsticepy.FieldPF
   :members:
   :undoc-members:

Find and run Solstice programs
-------------------------------------

.. autofunction:: solsticepy.find_solstice_root
.. autofunction:: solsticepy.find_prog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
