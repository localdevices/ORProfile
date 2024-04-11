.. _install:

============
Installation
============

After you have done a survey you should have obtained at least one, and (with a drone) two datasets. ORProfile offers
an Application Programming Interface (API) for loading the surveyed points, loading a point cloud derived through the
photogrammetric survey with drone, and combining these into one bathymetric dataset. At the moment we don't have a user
interface yet, but the API  is generic enough to establish this in a later stage, or offer this in an existing
environment such as QGIS.

Prerequisites
-------------

- a laptop or desktop computer with windows or linux
- python installed (version 3.9 or above)
- a cup of coffee or tea for awaiting finalization of the installtion

Installation works best with pip. We recommend to make a virtual environment first before installing the software.
The full set of commands for linux is:

.. code-block:: shell

    python -m venv $HOME/venv/orprofile
    source ${HOME}/venv/orprofile/bin/activate
    pip install orprofile

in windows you can do a similar thing:

.. code-block:: shell

    python -m venv %HOME%/venv/orprofile
    %HOME%/venv/orprofile/bin/activate
    pip install orprofile

If you wish to use the notebooks as example, then please grab a copy of our source code and look in the notebooks
folder.
