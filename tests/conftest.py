import os
from pytest import fixture

import geopandas as gpd

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
def splines(fn_splines):
    # read geojson in memory and return
    gdf = gpd.read_file(fn_splines)
    return gdf

# @fixture
# def mesh(splines):
