import matplotlib.pyplot as plt
from orprofile.api import depth_sample
def test_regrid_samples(mesh_points):

    ds = depth_sample.regrid_samples(mesh_points, mesh_points.points)
    # f = plt.figure(figsize=(16, 9))
    ax = mesh_points.plot()
    ds = ds.ugrid.to_crs(ax.projection)
    ds["samples"].ugrid.plot(ax=ax, alpha=0.8, zorder=4) #, alpha=0.8)
    plt.show()
