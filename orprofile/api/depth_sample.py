import numpy as np

def regrid_samples(mesh, gdf, nodata=-9999, column="Depth"):
    """
    regrid a set of depth samples to the provided mesh object

    Parameters
    ----------
    mesh

    Returns
    -------

    """
    # get the index numbers of each mesh cell
    da_index = mesh._get_index_da()
    # get index number of each sample
    points = da_index.ugrid.sel_points(
        x=gdf.geometry.x,
        y=gdf.geometry.y
    )
    face_index = np.ones(len(gdf), dtype=np.int64) * nodata  # -9999 is the nodata value
    face_index[points.mesh2d_index] = points.values
    gdf["face_index"] = face_index

    depth_mean = gdf[gdf["face_index"] != nodata][[column, "face_index"]].groupby("face_index").mean()
    depth_grid = np.zeros(da_index.shape) * np.nan
    depth_grid[depth_mean.index] = depth_mean[column]
    ds = mesh._get_empty_ds()
    ds["samples"] = ("mesh2d_nFaces", depth_grid)
    return ds

