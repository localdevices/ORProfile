
================================
OpenRiverProfile - Documentation
================================

OpenRiverProfile offers hardware setups and software methods to measure very precise bathymetry with only a smartphone,
a fishfinder and a low-cost RTK GNSS receiver and antenna. With a simple drone, also dry banks can be surveyed and added
to your project. As a result you will get a curvilinear mesh or GeoTIFF terrain model of the river or

For alluvial rivers, we provide means to estimate the total conveyance capacity of the stream and shape parameters of
the stream using typical geomorphological principles.

For small reservoirs and lakes, we offer geospatial interpolation to fill in gaps and provide a seamless bathymetric
model.

At present the software can be operated through the use of a jupyter notebook. We recommend to use jupyter instead of
a fixed script because some parameters can be refined using interactive plots.

This software is AGPL-3 licensed. Please read our LICENSE file for further information.

.. note::

    Acknowledgement: this software has been created within the OPEN PROFILE project, funded by the World Meteorological
    Organisation.

FIGURE!!

:ref:`Hardware design <hardware>`
:ref:`Guideline for bathymetric survey <survey>`
:ref:`Installation software <install>`
:ref:`Example notebook for bathymetry <notebook>`

.. toctree::
   :titlesonly:
   :hidden:
   :maxdepth: 2

   hardware.rst
   survey.rst
   install.rst
   notebook.rst
   api.rst
