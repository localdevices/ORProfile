from orprofile.api import Mesh
import matplotlib.pyplot as plt
import matplotlib


def test_mesh(splines):
    # setup a fresh mesh
    matplotlib.use("Qt5Agg")
    mesh = Mesh(splines, n=20, m=20)
    print(mesh)
    mesh.plot()
    plt.show()


def test_mesh_points(splines, fn_points):
    # setup a fresh mesh
    matplotlib.use("Qt5Agg")
    mesh = Mesh(splines, n=20, m=20, points=fn_points)
    print(mesh)
    mesh.plot()
    plt.show()
