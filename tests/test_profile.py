import orprofile

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def test_get_dist(mesh):
    da = mesh.map_rowwise(
        orprofile.profile.get_dist,
        name="distance from bank"
    )
    da.ugrid.plot()
    plt.show()


def test_get_depth(mesh_high_res):
    da = mesh_high_res.map_rowwise(
        orprofile.profile.get_depth,
        name="depth from bank",
        phi_left=0.05*np.pi,
        phi_right=0.25*np.pi
    )
    da.ugrid.plot()
    plt.show()
