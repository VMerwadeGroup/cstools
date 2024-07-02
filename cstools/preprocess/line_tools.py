from __future__ import annotations

import numpy as np
from scipy import interpolate
from shapely.geometry import LineString, Point, MultiPoint, Polygon
from shapely.ops import voronoi_diagram, nearest_points
import geopandas as gpd

def create_cl_from_bd(bd_geom: Polygon) -> LineString:

    ## TODO: 
    ## The algorithm will have troubles some shapes.
    ## There will be many Y-branches in the voronoi_diagram. 

    v_diag = voronoi_diagram(bd_geom, tolerance=0)
    v_geom = [poly for poly in v_diag.geoms]
    v_gdf = gpd.GeoDataFrame(geometry=v_geom)
    v_clip_gdf = gpd.clip(v_gdf, bd_geom).boundary

    ## create adjacent table
    bd_geom_small = bd_geom.buffer(-1e-2)
    points = {}
    for line in v_clip_gdf.geometry:
        coords = [c for c in line.coords]
        for c1, c2 in zip(coords[:-1], coords[1:]):
            p1, p2 = Point(c1), Point(c2)
            if (p1.within(bd_geom_small)):
                if (p1 not in points):
                    points[p1] = set()
                if (p2.within(bd_geom_small)):
                    points[p1].add(p2)
   

    ## check if there is any branches in the graph
    end_points = []
    split_points = []
    for point, adj in points.items():
        if (len(adj) == 1):     end_points.append(point)    
        elif (len(adj) > 2):    split_points.append(point)

    if (len(end_points) == 2):
        split_points = end_points
        
    current = split_points[0]
    visited = {current}
    path = [current]
    while (current != split_points[-1]):
        adjs = points[current]
        for adj in adjs:
            if adj not in visited:
                next = adj
                current = next
                visited.add(current)
                path.append(current)
                break

        if (len(end_points) > 2) and (current in end_points):
            current = split_points[0]
            path = [current]

    if (len(end_points) > 2):
        end_points1 = nearest_points(end_points[0], MultiPoint(end_points[1:]))
        end_points2 = [end_point for end_point in end_points if end_point not in end_points1]

        def create_middle_point(points):
            p1, p2 = points
            x1, y1 = p1.xy
            x2, y2 = p2.xy
            return Point(np.mean(x1+x2), np.mean(y1+y2))

        end1 = create_middle_point(end_points1)
        end2 = create_middle_point(end_points2)
        path = [end1] + path + [end2]   

    return LineString(path)


def smooth_line_geom(line_geom: LineString, spacing: float|None = 10) -> LineString:
    num = int(line_geom.length / spacing) + 1
    xy = np.array(line_geom.coords)
    xy = xy[:, :2]
    ## Note: splprep and splrep use different way to get the default s 
    tck, u = interpolate.splprep(xy.T, u=None, k=3, s=num) 
    u_new = np.linspace(u.min(), u.max(), num)
    new_x, new_y = interpolate.splev(u_new, tck)
    return LineString(zip(new_x, new_y))


    
