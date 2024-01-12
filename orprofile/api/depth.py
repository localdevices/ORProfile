# DataArray functions for estimating 2D depth

import numpy as np
import xarray as xr

from ..profile.lacey import depth_profile_left_right_parameters, deepest_point_y_dist
from .mesh import get_dist

def _get_gamma(xi, yi, c, alpha, L, kappa, A, h_m):
    """
    Estimate the asymmetry parameter for an entire curvilinear grid using the bank-to-bank distances (in yi)
    and the longitudinal distances (in xi).

    Parameters
    ----------
    xi : xugrid DataArray (float, holding rows/cols properties)
        longitudinal distances along river stretch [m]
    yi : array-like  (float, holding rows/cols properties)
        perpendicular distances (bank-to-bank) [m]
    c : float
        levee-to-levee / grid width ratio [-]
    alpha : float
        amplitude of the offset from the center [m]
    L : float
        wave length of meander shape [m]
    kappa : float [-np.pi, np.pi],
        phase offset of meander shape [rad]
    A : float
        conveyance, measured as wetted cross sectional surface [m2]
    h_m : array-like (float)
        depth at deepest point

    Returns
    -------
    xugrid DataArray (float, holding rows/cols properties)
        gamma parameter, describing the channel asymmetry
    """
    # compute the location of the middle ordinate coordinate
    y_mean = yi.groupby("cols").mean()
    y_disti = deepest_point_y_dist(xi, alpha, L, kappa, y_off_mean=0.)  # y_off_mean could possibly become an extra parameter

    h = 2 * h_m / np.pi
    Br_vector = y_mean - yi.groupby("cols").min() + y_disti.groupby("cols").max()
    Br_grid = np.zeros(len(xi))
    for n, br in enumerate(Br_vector.values):
        Br_grid[xi["cols"] == n] = br

    # # compute the grid width (B tilde)
    # B_grid = yi[-1] - yi[0]
    # Br_grid = B_grid - (y_mean - yi[0] + y_disti.max(axis=0))

    gamma = A/(h*Br_grid*c) - 1
    # gamma = y_disti/(c*B_grid) + 1
    # then we compute how far in the transect we are per grid cell
    return gamma

def _depth_profile_row(row, A):
    """
    solves Lacey equations with asymmetry over one row with already calculated gamma, h_m

    Parameters
    ----------
    row : xr.Dataset
        one single row of a full curvilinear grid, holding all information to calculate depth
    A : float
        cross-sectional surface [m2]
    c : float
        levee-to-levee / grid width ratio [-]

    Returns
    -------
    depth : xr.DataArray
        depth over the selected row of data, with the derived parameters
    """
    s = row.yi
    gamma = row["gamma"]
    h_m = row["h_m"]
    depth = depth_profile_left_right_parameters(s, gamma, h_m, A) * (h_m * 0 + 1)
    return depth


def depth_2d(mesh, alpha, L, kappa, c, A):
    """
    Lacey, solved 2-dimensionally, wrapped in xugrid data model

    Parameters
    ----------
    mesh : api.Mesh
        curvilinear mesh, to use to apply Lacey on.
    c : float
        levee-to-levee / grid width ratio [-]
    A : float
        cross-sectional surface [m2]
    alpha : float
        amplitude of the offset from the center [m]
    L : float
        wave length of meander shape [m]
    kappa : float [-np.pi, np.pi],
        phase offset of meander shape [rad]

    Returns
    -------
    """
    ds = mesh._get_empty_ds()
    xi = mesh.map_rowcol_wise(
        get_dist,
        name="xi",
        rowcol="rows"
    )

    yi = mesh.map_rowcol_wise(
        get_dist,
        name="yi",
        rowcol="cols"
    )
    ds["yi"] = yi
    ds["xi"] = yi

    # get the asymmetry parameter
    B_vector = ds["yi"].groupby("cols").max() - ds["yi"].groupby("cols").min()
    B_grid = np.zeros(len(ds["xi"]))
    for n, b in enumerate(B_vector.values):
        B_grid[ds["cols"] == n] = b

    # gamma = get_gamma(xi, yi, c, alpha, L, kappa)
    # get the h_m values
    h_m = A * np.pi / (2 * c * B_grid)
    gamma = _get_gamma(xi, yi, c, alpha, L, kappa, A, h_m)
    ds["h_m"] = ("mesh2d_nFaces", h_m)
    ds["gamma"] = ("mesh2d_nFaces", gamma)
    ds_g = xr.Dataset(ds).groupby("cols")
    depth = ds_g.map(
        _depth_profile_row,
        A=A,
    )
    depth.name = "depth"
    ds["depth"] = depth
    return ds
