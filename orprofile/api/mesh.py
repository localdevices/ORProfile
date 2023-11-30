import matplotlib.pyplot as plt
import numpy as np
import os

import geopandas as gpd

from meshkernel import (
    CurvilinearParameters,
    MakeGridParameters,
    GeometryList,
    MeshKernel,
    SplinesToCurvilinearParameters,
    OrthogonalizationParameters
)

class Mesh(object):
    def __init__(
            self,
            spline_shape: gpd.GeoDataFrame,
            n=10,
            m=10,
    ):
        self.n = n
        self.m = m
        self.splines = spline_shape

    @property
    def splines(self):
        return self._splines

    @splines.setter
    def splines(self, spline_shape):
        # dependent on what spline_shape is, either set it directly or read from file
        if isinstance(spline_shape, gpd.GeoDataFrame):
            self._splines = spline_shape
        else:
            if os.path.isfile(spline_shape):
                self._splines = gpd.read_file(spline_shape)


    @property
    def _is_splines_crossing(self):
        # TODO: check if each spline is crossing another line
        # if ...
        return True
        # else:
        #     return False
        # raise NotImplementedError

    @property
    def splines_mesh(self):
        if self._is_splines_crossing:
            separator = -999.0
            # loop over geometry and make a meshkernel compatible GeometryList object
            splines_x = []
            splines_y = []
            for n, geom in enumerate(self.splines.geometry):
                if n == 0:
                    splines_x += list(geom.xy[0])
                    splines_y += list(geom.xy[1])
                else:
                    splines_x += [separator] + list(geom.xy[0])
                    splines_y += [separator] + list(geom.xy[0])

            # splines_x = np.array([2.0, 4.0, 7.0, separator,
            #                       -1.0, 1.0, 5.0, separator,
            #                       3.0, -2.0, separator,
            #                       7.0, 4.0], dtype=np.double)
            # splines_y = np.array([1.0, 3.0, 4.0, separator,
            #                       4.0, 6.0, 7.0, separator,
            #                       1.0, 6.0, separator,
            #                       3.0, 8.0], dtype=np.double)
            splines_x = np.float64(splines_x)
            splines_y = np.float64(splines_y)
            return GeometryList(splines_x, splines_y)
        else:
            print("Splines are not crossing each other properly, ensure each spline crosses at least 2 other splines")

        raise NotImplementedError
        def plot(ax=None):
            raise NotImplementedError

    @property
    def mesh_kernel(self):
        mk = MeshKernel()
        curvilinear_parameters = CurvilinearParameters()
        curvilinear_parameters.n_refinement = self.n
        curvilinear_parameters.m_refinement = self.m
        mk.curvilinear_compute_transfinite_from_splines(
            self.splines_mesh, curvilinear_parameters
        )
        return mk.curvilineargrid_get()

    def plot(self, ax=None):
        if ax is None:
            f, ax = plt.subplots()
        self.mesh_kernel.plot_edges(ax)

