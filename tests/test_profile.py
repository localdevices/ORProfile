import orprofile

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import matplotlib

def test_get_dist(mesh):
    da = mesh.map_rowcol_wise(
        orprofile.profile.get_dist,
        name="y",
        rowcol="cols"
    )
    da.ugrid.plot()
    plt.show()


def test_get_depth(mesh_high_res):
    ax = mesh_high_res.plot()
    da = mesh_high_res.map_rowcol_wise(
        orprofile.profile.get_depth,
        name="depth from bank",
        phi_left=0.15*np.pi,
        phi_right=0.005*np.pi
    )
    da_proj = da.ugrid.to_crs(ax.projection)
    da_proj.ugrid.plot(ax=ax, vmax=0., alpha=0.4, cmap="viridis", zorder=3)  #ax=ax, transform=ccrs.epsg(mesh_high_res.splines.crs.to_epsg()),
    plt.show()
