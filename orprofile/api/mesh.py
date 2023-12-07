import matplotlib.pyplot as plt
import numpy as np
import os

import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import geopandas as gpd
from shapely.geometry import Point

from meshkernel import (
    CurvilinearParameters,
    MakeGridParameters,
    GeometryList,
    MeshKernel,
    SplinesToCurvilinearParameters,
    OrthogonalizationParameters
)
from xugrid import Ugrid2d

class Mesh(object):
    def __init__(
            self,
            spline_shape: gpd.GeoDataFrame,
            n=10,
            m=10,
            points=None
    ):
        self.n = n
        self.m = m
        self.splines = spline_shape
        if points is not None:
            self.points = points

    def __repr__(self):
        info = {
            "n": self.n,
            "m": self.m,
            "crs": self.crs,
            "points head": None if self.points is None else self.points["geometry"].head()
        }
        info_str = ""
        for k, v in info.items():
            info_str += f"{k}: {v}\n"
        return info_str


    @property
    def crs(self):
        return self.splines.crs

    @property
    def points(self):
        if hasattr(self, "_points"):
            return self._points

    @points.setter
    def points(self, gdf):
        if isinstance(gdf, gpd.GeoDataFrame):
            self._points = gdf
        else:
            self.read_points(gdf)

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
                    splines_y += [separator] + list(geom.xy[1])
            # convert into 64-bit floats (needed for kernel)
            splines_x = np.float64(splines_x)
            splines_y = np.float64(splines_y)
            return GeometryList(splines_x, splines_y)
        else:
            print("Splines are not crossing each other properly, ensure each spline crosses at least 2 other splines")

        raise NotImplementedError
    def plot(
            self,
            ax=None,
            tiles="GoogleTiles",
            extent=None,
            zoom_level=18,
            tiles_kwargs={"style": "satellite"},
            **kwargs
    ):
        if tiles is not None:
            tiler = getattr(cimgt, tiles)(**tiles_kwargs)
            crs = tiler.crs
        else:
            crs = ccrs.PlateCarree()
        if ax is None:
            ax = plt.subplot(projection=crs)
        if extent is not None:
            ax.set_extent(extent, crs=ccrs.PlateCarree())
        m = self.mesh2d.to_crs(crs.to_wkt())
        m.plot(ax=ax, label="mesh", zorder=2)
        if tiles is not None:
            ax.add_image(tiler, zoom_level, zorder=1)
        # now add the gdf
        # self.mesh_kernel.plot_edges(ax, color="r", transform=ccrs.epsg(self.splines.crs.to_epsg()), label="mesh edges")
        self.splines.plot(ax=ax, transform=ccrs.epsg(self.splines.crs.to_epsg()), zorder=2, color="c", linewidth=2., label="splines")
        if self.points is not None:
            self.points.plot(
                column="depth",
                ax=ax,
                transform=ccrs.epsg(self.splines.crs.to_epsg()),
                zorder=2,
                # color=self.points.geometry.z,
                marker="+",
                markersize=30,
                # edgecolor="w",
                label="sonar survey",
                legend=True
            )

        # ax.add_geometries(gdf["geometry"], crs=ccrs.PlateCarree())
        ax.legend()
        return ax

    @property
    def mesh_kernel(self):
        mk = MeshKernel()
        curvilinear_parameters = CurvilinearParameters()
        curvilinear_parameters.n_refinement = self.n
        curvilinear_parameters.m_refinement = self.m
        mk.curvilinear_compute_transfinite_from_splines(
            self.splines_mesh, curvilinear_parameters
        )
        mk.curvilinear_convert_to_mesh2d()  #.curvilineargrid_get()
        return mk


    @property
    def mesh2d(self):
        grid = self.mesh_kernel.mesh2d_get()
        return Ugrid2d.from_meshkernel(grid, crs=self.crs)



    def read_points(self, fn, crs=None):
        """
        Read a set of surveyed points. These should be contained in the point geometry of the read file

        Parameters
        ----------
        fn

        Returns
        -------

        """
        gdf = gpd.read_file(fn)
        # check if all geometries are points
        assert(
            all(
                [isinstance(geom, Point) for geom in gdf.geometry]
            )
        ), f'shapefile may only contain geometries of type "Point"'
        if not(hasattr(gdf, "crs")):
            if crs is None:
                raise ValueError(
                    "a CRS must be provided either within the provided point geomtry file or explicitly "
                    "through `crs=`"
                )
            gdf.set_crs(crs)
        gdf["depth"] = gdf.geometry.z
        self.points = gdf