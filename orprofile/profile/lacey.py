"""contains functions that are mapped over profiles (i.e. cross sections, perpendicular to flow directions"""
import numpy as np


def deepest_point_y_dist(x, alpha, L, kappa, y_off_mean=0.):
    """
    Estimates the location of the deepest point as offset from the centerline
    as function of one or more longitudinal
    coordinates x, assuming a sinusoidal function. If x is in meters, then the result
    is also in meters.

    Parameters
    ----------
    x : float, array
        one or multiple x coordinates (longitudinal direction). can also be a grid
        of coordinates [m]
    alpha : float,
        amplitude of the offset from the center [m]
    L : float
        wave length of meander shape [m]
    kappa : float [-np.pi, np.pi],
        phase offset of meander shape [rad]
    y_off_mean : float, optional
        centerline offset

    Returns
    -------
    float, array
        distance of deepest point from centerline
    """
    return y_off_mean + alpha * np.sin(2 * np.pi * (x / L) - kappa)


def depth_from_width(s, phi=0.25 * np.pi, h_m=10):
    """
    Basic implementation of Eq. 12 from Savenije 2003, where the first coordinate (zero meters) is the lowest point.
    We cut off numbers at zero as that is where we reach the natural levee crest.

    Parameters
    ----------
    s : float, vector-like
        ordinate coordinates, measured as distance from the deepest point [m]
    phi : float
        natural angle of repose of the sediment on the natural levees [rad]
    h_m : float
        depth at deepest point

    Returns
    -------
    float, vector-like
        depth (negative numbers only)
    """
    # zero is the lowest point in the stream starting in the middle of the stream
    angle = np.minimum(np.tan(phi) / h_m * s, 0.5 * np.pi)
    depth = -np.cos(angle) * h_m
    return depth


def depth_profile_left_right_parameters(s, gamma, h_m, A):
    """
    Savenije eq. 12 applied with asymmetrical channel location, defined by parameter gamma.

    The implementation assumes that the most extreme coordinates are at more or less the same elevation and
    within the area spanning between both natural levees (they do not have to be at the natural levee per se)

    Parameters
    ----------
    s : float, vector-like
        cross-section coordinates measured from zero (left-bank-side) to maximum (right-bank-side)
    gamma : float
        channel asymmetry coefficient [-]
    h_m : float
        depth at deepest point
    A : float
        conveyance, measured as wetted cross sectional surface [m2]

    Returns
    -------
    float, vector-like
        depth over left to right bank side coordinates (negative numbers only)

    """
    B = s[-1] - s[0]
    frac_right = 1 / (1 + gamma)

    # make coordinates with zero at the deepest point negative (positive) left (right) of deepest point
    s_ordinal = s - (1 - frac_right) * B
    # estimate the average depth (eq. 15 in Savenije 2003)
    h = 2 * h_m / np.pi
    # solve B_r and B_l
    B_r = A / (h * (1 + gamma))
    B_l = gamma * B_r
    # compute tan phi left and right
    tan_phi_left = h_m * np.pi / (2 * B_l)
    tan_phi_right = h_m * np.pi / (2 * B_r)
    tan_phi_right2 = gamma * tan_phi_left
    phi_left = np.arctan(tan_phi_left)
    phi_right = np.arctan(tan_phi_right)

    depth = np.zeros(s.shape)

    depth[s_ordinal < 0] = depth_from_width(-s_ordinal[s_ordinal < 0], phi=phi_left[s_ordinal < 0],
                                            h_m=h_m[s_ordinal < 0])
    depth[s_ordinal >= 0] = depth_from_width(s_ordinal[s_ordinal >= 0], phi=phi_right[s_ordinal >= 0],
                                             h_m=h_m[s_ordinal >= 0])
    return np.minimum(depth, 0)
