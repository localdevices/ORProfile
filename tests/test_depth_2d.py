import numpy as np
import orprofile
import matplotlib.pyplot as plt

def test_depth(mesh_high_res):
    alpha = 20
    L = 150
    kappa = 0.5*np.pi
    c = 0.8
    A = 500
    ds = orprofile.api.depth_2d(mesh_high_res, alpha, L, kappa, c, A)
    ax = plt.subplot(111)
    p = ds["depth"].ugrid.plot(ax=ax)
    plt.show()
