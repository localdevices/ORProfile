from orprofile.api import Mesh
import matplotlib.pyplot as plt

def test_mesh(splines):
    # setup a fresh mesh
    mesh = Mesh(splines)
    print(mesh)
    mesh.plot()
    plt.show()