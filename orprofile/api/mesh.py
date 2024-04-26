import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import copy
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import os
from shapely.geometry import Point
from meshkernel import (
    CurvilinearParameters,
    GeometryList,
    MeshKernel,
)
import xarray as xr
import xugrid as xu

__all__ = ["Mesh"]

def get_dist(row):
    # get the centroid coordinates of each face
    x = row.ugrid.to_geodataframe().centroid.x
    y = row.ugrid.to_geodataframe().centroid.y
    # compute the difference in distance per grid cell from one bank to the other
    ds = (np.diff(x)**2 + np.diff(y)**2)**0.5
    # distance from bank
    s = np.cumsum(np.pad(ds, (1, 0), "constant"))
    # use the "cols" variable as a template
    return xr.DataArray(s, coords=row.coords)

def get_xi(mesh):
    return mesh.map_rowcol_wise(
        get_dist,
        name="xi",
        rowcol="rows"
    )

def get_yi(mesh):
    return mesh.map_rowcol_wise(
        get_dist,
        name="yi",
        rowcol="cols"
    )


def map_func(row, f, ugrid, name="new", **kwargs):
    grid_sel = ugrid.isel(mesh2d_nFaces=row.mesh2d_nFaces.values)
    row_ug = xu.UgridDataset(row, grids=grid_sel)
    # now apply the function
    result = f(row_ug, **kwargs)
    result.name = name
    return result

# specific function that computes something. This can be made by the user and applied

class Mesh(object):
    def __init__(
            self,
            spline_shape: gpd.GeoDataFrame,
            n=10,
            m=10,
            points=None,
            mesh_kernel=None  # if provided it will be set from mesh kernel, otherwise grown from the pars and splines
    ):
        self.n = n
        self.m = m
        self.splines = spline_shape
        if points is not None:
            self.points = points
        self.mesh_kernel = mesh_kernel


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
        """

        Returns
        -------
        gpd.GeoDataFrame
            loaded sonar survey points
        """
        if hasattr(self, "_points"):
            return self._points

    @points.setter
    def points(self, gdf):
        if isinstance(gdf, gpd.GeoDataFrame):
            self._points = gdf
        else:
            self.read_points(gdf)


    @property
    def mesh_kernel(self):
        """

        Returns
        -------
        meshkernel.MeshKernel
            contains the curvilinear grid based on defined splines and row/column parameters
        """
        if hasattr(self, "_mesh_kernel"):
            return self._mesh_kernel

    @mesh_kernel.setter
    def mesh_kernel(self, mesh_kernel=None):
        if not mesh_kernel:
            mesh_kernel = MeshKernel()
            curvilinear_parameters = CurvilinearParameters()
            curvilinear_parameters.n_refinement = self.n
            curvilinear_parameters.m_refinement = self.m
            mesh_kernel.curvilinear_compute_transfinite_from_splines(
                self.splines_mesh, curvilinear_parameters
            )
            mesh_kernel.curvilinear_convert_to_mesh2d()  #.curvilineargrid_get()
        self._mesh_kernel = mesh_kernel


    @property
    def splines(self):
        """

        Returns
        -------
        gpd.GeoDataFrame
            4 splines defining the area of interest (2 along banks and 2 perpendicular to river
        """
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
        """

        Returns
        -------
        A GeometryList that can be used to construct a mesh in meshkernel
        """
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



    @property
    def mesh2d(self):
        """

        Returns
        -------
        xu compatible 2D grid
        """
        grid = self.mesh_kernel.mesh2d_get()
        return xu.Ugrid2d.from_meshkernel(grid, crs=self.crs)


    @property
    def rowscols(self):
        rows = np.repeat([np.arange(self.n)], self.m, axis=0).flatten()
        columns = np.repeat(np.arange(self.m), self.n)
        return rows, columns

    def plot(
            self,
            ax=None,
            tiles="GoogleTiles",
            extent=None,
            zoom_level=18,
            tiles_kwargs={"style": "satellite"},
            splines_kw={},
            points_kw={},
            plot_points=True
    ):
        """
        Plot available data in a geographically aware plot (cartopy)

        Parameters
        ----------
        ax : plt.Axes, optional
            pre-defined axes (if required)
        tiles : str,
            cartopy.io.img_tiler submethod tiles set name to use as background (default: GoogleTiles)
        extent : list[float],
            extent as xmin, xmax, ymin, ymax to reduce the image to
        zoom_level : int, optional
            zoom level to use for tiles set (default: 18)
        tiles_kwargs : dict, optional
            kwargs to pass to cartopy.io.img_tiler.<tiles>. Default: {"style": "satellite"}
        kwargs

        Returns
        -------

        """
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
        m.plot(ax=ax, label="mesh", zorder=2, alpha=0.5)
        if tiles is not None:
            ax.add_image(tiler, zoom_level, zorder=1)
        # now add the gdf
        # self.mesh_kernel.plot_edges(ax, color="r", transform=ccrs.epsg(self.splines.crs.to_epsg()), label="mesh edges")
        self.splines.plot(
            ax=ax,
            transform=ccrs.epsg(self.splines.crs.to_epsg()),
            zorder=2,
            color="c",
            linewidth=2.,
            label="splines",
            **splines_kw
        )
        if self.points is not None and plot_points:
            self.points.plot(
                column="depth",
                ax=ax,
                transform=ccrs.epsg(self.splines.crs.to_epsg()),
                zorder=2,
                marker="+",
                markersize=30,
                label="sonar survey",
                legend=True,
                **points_kw
            )
        ax.legend()
        return ax


    def _get_empty_ds(self):
        """
        make an empty template for a xu.UgridDataset with a underlying unstructured grid

        Returns
        -------
        ds : xu.UgridDataset or xu.UgridDataArray object

        """

        # make a columns and row coordinate for each node
        rows = np.repeat([np.arange(self.n)], self.m, axis=0).flatten()
        columns = np.repeat(np.arange(self.m), self.n)

        # prepare UgridDataArrays, using the grid
        da_rows = xu.UgridDataArray(
            xr.DataArray(
                data=rows,
                dims=[self.mesh2d.face_dimension]
            ),
            grid=self.mesh2d,
        )
        da_cols = xu.UgridDataArray(
            xr.DataArray(
                data=columns,
                dims=[self.mesh2d.face_dimension]
            ),
            grid=self.mesh2d,
        )

        # combine everything into a UgridDataSet
        ds = xu.UgridDataset(grids=self.mesh2d)
        ds.coords["rows"] = da_rows
        ds.coords["cols"] = da_cols
        ds["rows"] = da_rows
        ds["cols"] = da_cols


        return ds


    def _get_index_da(self):
        """
        Retrieve a xu.UgridDataArray with the index as data

        Returns
        -------

        """

        return xu.UgridDataArray(
            xr.DataArray(
                data=self.mesh2d.to_dataset()["mesh2d_nFaces"],
                dims=[self.mesh2d.face_dimension]
            ),
            grid=self.mesh2d,
        )

    def map_rowcol_wise(self, func, name="new_var", rowcol="rows", **kwargs):
        """
        Map a 1D function over each row in the defined mesh object
        under self.mesh2d. Output is returned with the mesh intelligence
        as xu.UgridDataArray object.

        Parameters
        ----------
        func : def
            function to map over
        name : str, optional
            name of output xr.DataArray object
        kwargs
            : keyword arguments to pass to func

        Returns
        -------
        da : xu.UgridDataArray object

        """

        # create an empty ds with grid properties
        ds = self._get_empty_ds()
        ds_g = ds.groupby(rowcol)
        da = ds_g.map(
            map_func,
            f=func,
            name=name,
            ugrid=self.mesh2d,
            **kwargs
        )
        return xu.UgridDataArray(
            da,
            grid=self.mesh2d
        )


    def read_points(self, fn, crs=None):
        """
        Read a set of surveyed points. These should be contained in the point geometry of the read file

        Parameters
        ----------
        fn : str, path-like,
            geopandas compatible geographic vector file (e.g. GeoJSON)
        crs : str, int, optional
            coordinate reference system, can be anything that can be interpreted by pyproj.CRS.from_user_input

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


def mesh_new_row(mesh, side="left"):
    mesh_new = copy.copy(mesh)
    mesh2d = mesh.mesh2d.to_dataset()
    rows, _ = mesh.rowscols
    if side == "left":
        bank_idx = rows == rows.min()
        first_extend_idx = [3, 0]
        extend_idx = [2, 1]
    else:
        bank_idx = rows == rows.max()
        first_extend_idx = [1, 2]
        extend_idx = [0, 3]
    mesh2d_faces_row = mesh2d["mesh2d_face_nodes"][bank_idx]
    new_x, new_y = [], []
    new_edge = []
    cur_nr_nodes = len(mesh2d["mesh2d_node_x"])
    for n, cell in enumerate(mesh2d_faces_row):
        x_sel = mesh2d["mesh2d_node_x"][cell].values
        y_sel = mesh2d["mesh2d_node_y"][cell].values
        idx_sel = mesh2d["mesh2d_nNodes"][cell].values

        if n == 0:
            # select the nodes that lie on the edge of the domain
            x1, x2 = x_sel[first_extend_idx]
            y1, y2 = y_sel[first_extend_idx]
            _, idx_edge1 = idx_sel[first_extend_idx]
            new_x.append(x2 - (x1 - x2))
            new_y.append(y2 - (y1 - y2))
            new_edge.append([idx_edge1, cur_nr_nodes])
            cur_nr_nodes += 1
        x1, x2 = x_sel[extend_idx]
        y1, y2 = y_sel[extend_idx]
        _, idx_edge2 = idx_sel[extend_idx]
        new_x.append(x2 - (x1 - x2))
        new_y.append(y2 - (y1 - y2))
        new_edge.append([idx_edge2, cur_nr_nodes])
        # connect hanging nodes
        new_edge.append([cur_nr_nodes - 1, cur_nr_nodes])
        # make the last edge the first for the next cell
        idx_edge1 = copy.deepcopy(idx_edge2)
        cur_nr_nodes += 1
    # now generate the new nodes and faces
    x = np.concatenate([mesh2d["mesh2d_node_x"].values, new_x])
    y = np.concatenate([mesh2d["mesh2d_node_y"].values, new_y])
    new_faces = np.vstack(
        [mesh2d_faces_row.values[:, 0],
         mesh2d_faces_row.values[:, 1],
         np.arange(len(mesh2d_faces_row)) + len(mesh2d["mesh2d_node_x"]),
         np.arange(len(mesh2d_faces_row)) + 1 + len(mesh2d["mesh2d_node_x"])]
    ).T
    faces = np.insert(mesh2d["mesh2d_face_nodes"], np.where(bank_idx)[0], new_faces, axis=0)
    print(len(faces))
    edges = np.vstack([mesh2d["mesh2d_edge_nodes"], new_edge])
    # now rearrange the mesh
    del mesh2d["mesh2d_node_x"]
    del mesh2d["mesh2d_node_y"]
    mesh2d["mesh2d_nNodes"] = np.arange(len(x), dtype=np.int64)
    mesh2d["mesh2d_node_x"] = "mesh2d_nNodes", x
    mesh2d["mesh2d_node_y"] = "mesh2d_nNodes", y
    mesh2d["mesh2d_node_x"].attrs = mesh.mesh2d.to_dataset()["mesh2d_node_x"].attrs
    mesh2d["mesh2d_node_y"].attrs = mesh.mesh2d.to_dataset()["mesh2d_node_y"].attrs
    # first remove the original faces
    del mesh2d["mesh2d_face_nodes"]
    # then replace the number of faces with the new length
    mesh2d["mesh2d_nFaces"] = np.arange(len(faces), dtype=np.int64)
    # now reinsert the new list of faces
    mesh2d["mesh2d_face_nodes"] = faces
    # do the same for the edges
    del mesh2d["mesh2d_edge_nodes"]
    mesh2d["mesh2d_nEdges"] = np.arange(len(edges), dtype=np.int64)
    mesh2d["mesh2d_edge_nodes"] = ("mesh2d_nEdges", "two"), edges
    new_mesh = xu.Ugrid2d.from_dataset(mesh2d)
    mesh_new.mesh_kernel = new_mesh.meshkernel
    # increase amount of rows
    mesh_new.n += 1
    # print(mesh2d)
    return mesh_new
