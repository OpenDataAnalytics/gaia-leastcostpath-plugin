Gaia Least Cost Path Plugin
================================

This is a plugin for Gaia (https://github.com/OpenDataAnalytics/gaia) that
calculates the least cost path between two points over a raster surface.
Least cost path analysis calculates the most cost-effective route between a source and destination.
Cost can be a function of elevation, time, or any other criteria that is represented as values on a raster grid, where a higher value indicates a higher cost.
As part of the analysis, the 8 neighbors of a grid cell are evaluated and the path moves to the cell with the smallest value.
This evaluation is repeated until the source and destination are connected.  The output is a vector line that connects the source and destination points.

An example of how to use this plugin can be found `here <gaia_processes.html>`__.

Installation
-----------------

- git clone https://github.com/OpenDataAnalytics/gaia-leastcostpath-plugin.git
- cd gaia-leastcostpath-plugin
- pip install -e .
- pip install -r requirements


Testing
-----------------

- pip install -r requirements-dev.txt
- python -m unittest discover


Expected Inputs
------------------
NOTE: This will change in the near future to accept standard GaiaIO class instances as inputs.

Currently the input must be a dict in the form of:

::

   {
      "uri": <filepath of raster image>,
      "start": (longitude, latitude),
      "end": (longitude, latitude)
   }



Table of Contents
-----------------

.. toctree::
   :maxdepth: 2

   gaia_leastcostpath


.. _Gaia: http://www.github.com/opendataanalytics/gaia
.. _Kitware: http://www.kitware.com
.. _Epidemico: http://epidemico.com
