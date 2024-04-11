.. currentmodule:: orprofile

.. _api:

=============
API reference
=============

The API of ``orprofile`` consists of several classes and functions. Here we document these in a structured way so that you
can reuse the code, build your own scripts, notebooks or even software around it.

.. note::

    We highly recommend to use the excellent xarray_ manual side-by-side with ``orprofile`` to understand
    how to effectively work with the xarray_ ecosystem.

In the remaining sections, we describe the API classes, and the functions they are based on.

.. _mesh:

Mesh class
==================


.. autosummary::
    :toctree: _generated

    api.Mesh
    api.Mesh.points
    api.Mesh.splines
    api.Mesh.mesh_kernel
    api.Mesh.mesh2d
    api.Mesh.plot
    api.Mesh.map_rowcol_wise
    api.Mesh.read_points


.. _depth:

Depth functions
===============

.. autosummary::
    :toctree: _generated


    api.depth_2d
    profile.depth_profile_left_right_parameters
    profile.depth_from_width
    profile.deepest_point_y_dist



.. _xarray: https://docs.xarray.dev/
.. _pandas: https://pandas.pydata.org/
