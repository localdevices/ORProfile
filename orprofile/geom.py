import geopandas as gpd
from shapely.ops import split, polygonize


def split_line_on_line(line1, line2):
    """
    split ``line1`` at crossing of ``line2`` or return None when no crossing is found.

    Parameters
    ----------
    line1 : LineString
        Line under consideration for splitting
    line2 : LineString
        Line to test for splitting point

    Returns
    -------
    List[LineString]
       contains two LineStrings, splitted from ``line1``.

    """
    intersection = line1.intersection(line2)
    if not intersection.is_empty:
        if not intersection.geom_type == "Point":
            raise ValueError(f"intersection must be of type Point but is {intersection.geom_type}")
        # intersection found! split the line
        return split(line1, line2)


def get_intersection_part(line1, other_lines):
    """
    Extract intersecting part of line between two other lines.

    Parameters
    ----------
    line1 : LineString
        Line under consideration for splitting
    other_lines: List[LineString]
        List of lines that are tested for crossing. There must be two unique other lines that cross ``line1`` for the
        function to work.

    Returns
    -------
    LineString
        the part of ``line1`` that falls in between two crossing lines found in ``other_lines``

    """
    split1 = False
    for n, line2 in enumerate(other_lines):
        if line2 != line1:
            # as soon as split1 is True something else should be done
            if not split1:
                split_lines = split_line_on_line(line1, line2)
                if split_lines is not None:
                    split1 = True
                    # now the split_lines should contain two lines
            else:
                # check if line goes through either one of the split_lines segments
                for l in split_lines.geoms:
                    split_line_again = split_line_on_line(l, line2)
                    if split_line_again is not None:
                        split_line_again = split(l, line2)
                        # check which part of the line fits with the former part
                        for split_l in split_line_again.geoms:
                            touch = split_l.buffer(0.001).overlaps(
                                split_lines.geoms[0].buffer(0.001)
                            ) and split_l.buffer(0.001).overlaps(
                                split_lines.geoms[1].buffer(0.001)
                            )
                            if touch:
                                return split_l


def crossing_lines_to_polygon(lines):
    """
    Convert 4 crossing lines into single enclosing polygon.

    Parameters
    ----------
    lines : list[LineString]
       lines that cross in a quasi-rectangular shape

    Returns
    -------
    Polygon

    """
    lines_cut = [get_intersection_part(l, lines) for l in lines]
    return gpd.GeoDataFrame(geometry=list(polygonize(lines_cut)))
