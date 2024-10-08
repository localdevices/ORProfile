from orprofile.api import Mesh
import matplotlib.pyplot as plt


def test_mesh(splines):
    # setup a fresh mesh
    # matplotlib.use("Qt5Agg")
    mesh = Mesh(splines, n=20, m=20)
    print(mesh)
    mesh.plot()
    plt.show()


def test_mesh_outside(splines):
    # setup a fresh mesh
    # matplotlib.use("Qt5Agg")
    mesh = Mesh(splines, n=20, m=20)
    mesh_ex = mesh.add_rows(left=5, right=5)
    # mesh_ex2 = Mesh(splines, n=mesh_ex.n, m=mesh_ex.m, mesh_kernel=mesh_ex.mesh_kernel)
    print(mesh)
    mesh_ex.plot()
    plt.show()


def test_mesh_points(splines, fn_points):
    # setup a fresh mesh
    # matplotlib.use("Qt5Agg")
    mesh = Mesh(splines, n=20, m=20, points=fn_points)
    print(mesh)
    mesh.plot()
    plt.show()

