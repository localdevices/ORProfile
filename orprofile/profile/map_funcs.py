"""contains functions that are mapped over profiles (i.e. cross sections, perpendicular to flow directions"""
import numpy as np
import xarray as xr

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



def depth_from_width(s, phi=0.1*np.pi, h_m=10):
    """
    Simplest implementation of Eq. 12 from Savenije 2003, where the first coordinate (zero meters) is the lowest point. We stop at zero as that
    is where we reach the natural levee crest.

    Parameters
    ----------
    s : array-like
        distances from middle (deepest point) to bank [m]
    phi : float, optional
        angle of repose [-]
    h_m : float, optional
        depth in middle

    Returns
    -------
    array-like
        depth at points evaluated over s

    """
    # zero is the lowest point in the stream starting in the middle of the stream
    angle = np.minimum(np.tan(phi)/h_m * s, 0.5*np.pi)
    depth = -np.cos(angle)*h_m
    return depth

def depth_profile_from_parameters(s, phi=None, h_m=None):
    """
    Solving the entire profile using a symmetrical assumption and either phi or h_m defined

    Parameters
    ----------
    s : array-like
        distances from left bank [m]
    phi : float, optional
        angle of repose [-], if not given, h_m must be provided
    h_m : float, optional
        depth in middle [m], if not given, phi must be provided

    Returns
    -------
    array-like
        depth evaluated over s

    """
    if phi is None and h_m is None:
        raise ValueError("Either phi or h_m has to be supplied. You supplied nothing.")
    if phi is not None and h_m is not None:
        raise ValueError("Either phi or h_m has to be supplied. You supplied both. Select either one of the two")
    # compute the total width of the channel
    B = s[-1] - s[0]
    if phi is None:
        # h_m was supplied, so resolve for phi
        phi = np.arctan(np.pi*h_m/B)
    else:
        # phi was supplied, so resolve for h_m
        h_m = B*np.tan(phi)/np.pi
    # compute ordinate coordinates (i.e. from center line towards banks)
    _s = s - 0.5*B
    return depth_from_width(_s, phi=phi, h_m=h_m)

def depth_profile_left_right_parameters(s, phi_left, phi_right):
    """
    Solving entiure profile with different angles at left and right bank. The solutions are weighted according
    to the relative distance from both banks.

    Parameters
    ----------
    s : array-like
        distances from left bank [m]
    phi_left : float, optional
        angle of repose at left bank [-]
    phi_right : float, optional
        angle of repose at right bank [-]

    Returns
    -------
    array-like
        depth evaluated over s [m]

    """
    depths_left = depth_profile_from_parameters(s, phi=phi_left)
    # right bank
    depths_right = depth_profile_from_parameters(s, phi=phi_right)
    weights = np.linspace(0, 1, len(s))
    return (1 - weights) * depths_left + weights * depths_right

def get_depth(row, phi_left, phi_right):
    """
    Computes distances from bank from coordinates provided in xr.DataArray row, and estimates depth using
    Savenije et al. 2003 with possibility of using asymmetrical depths

    Parameters
    ----------
    row : xu.UgridDataArray
        containing only one cross section
    phi_left : float
        left bank angle of repose [-]
    phi_right : float
        right bank angle of repose [-]

    Returns
    -------
    array-like
        depths over evaluated coordinates [m]
    """
    s = get_dist(row)
    # with s derive the cross section
    depth = depth_profile_left_right_parameters(s, phi_left, phi_right)
    return depth