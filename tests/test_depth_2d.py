import numpy as np
import orprofile
import matplotlib.pyplot as plt
import pytest


@pytest.mark.parametrize(
    "mesh_",
    [
        # "mesh_high_res",
        "mesh_extend",
     ]
)
def test_depth(mesh_, request):
    mesh_ = request.getfixturevalue(mesh_)
    alpha = 20
    L = 150
    kappa = 0.5*np.pi
    c = 0.8
    A = 500
    ds = orprofile.api.depth_2d(mesh_, alpha, L, kappa, c, A)
    ds = ds.ugrid.to_crs(4326)
    p = ds["depth"].ugrid.plot()
    plt.show()
