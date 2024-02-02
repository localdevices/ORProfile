import os
from pytest import fixture

import geopandas as gpd
from orprofile import api

EXAMPLE_DATA_DIR = os.path.join(
    os.path.split(__file__)[0], "..", "examples", "data"
)

@fixture
def fn_splines():
    return os.path.join(
        EXAMPLE_DATA_DIR,
        "splines.geojson"
    )


@fixture
def fn_points():
    return os.path.join(
        EXAMPLE_DATA_DIR,
        "bamboi_survey.geojson"
    )

@fixture
def splines(fn_splines):
    # read geojson in memory and return
    gdf = gpd.read_file(fn_splines)
    # gdf = gpd.GeoDataFrame(pd.concat([gdf.iloc[0], gdf.iloc[3], gdf.iloc[1], gdf.iloc[2]], axis=1).T)
    return gdf

@fixture
def mesh(splines):
    return api.mesh.Mesh(splines, n=10, m=20)

@fixture
def mesh_high_res(splines):
    return api.mesh.Mesh(splines, n=50, m=100)

@fixture
def mesh_points(splines, fn_points):
    # setup a fresh mesh
    # matplotlib.use("Qt5Agg")
    return api.mesh.Mesh(splines, n=10, m=20, points=fn_points)




# @fixture
# def mesh(splines):
