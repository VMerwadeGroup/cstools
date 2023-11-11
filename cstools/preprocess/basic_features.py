from __future__ import annotations

import numpy as np
from shapely.geometry import Point, Polygon, LineString
from shapely import ops
import geopandas as gpd

def get_approx_channel_length(bd_geom: Polygon) -> float:
    '''
    assume the boudary can be approximated to a rectangle
    the quadratic equation solutions provide length and width
    -b +- sqrt(b^2 - 4ac) / 2a, the positive one is length	
    '''

    p = bd_geom.length # perimeter
    a = bd_geom.area # area
    return (p/2 + np.sqrt(p*p/4 - 4*a))/2

def get_approx_channel_width(bd_geom: Polygon) -> float:
    '''
    assume the boudary can be approximated to a rectangle
    the quadratic equation solutions provide length and width
    -b +- sqrt(b^2 - 4ac) / 2a, the negative one is width
    '''
    
    p = bd_geom.length # perimeter
    a = bd_geom.area # area
    return (p/2 - np.sqrt(p*p/4 - 4*a))/2

def get_avg_channel_width(bd_geom: Polygon, cl_geom: LineString) -> float:
    '''
    Assume the river channel in the boundary can be approximated as a rectangle.
    Thus, the average channel width is area / length.
    '''
    cl_in_bd = bd_geom.intersection(cl_geom)
    area = bd_geom.area
    length = cl_in_bd.length
    return area / length

def get_cl_from_multilines(
        cl_gdf: gpd.GeoDataFrame, 
        bd_gdf: gpd.GeoDataFrame, 
        rootDir: str | None = None
    ) -> gpd.GeoDataFrame:

    ## preprocessing and clip line features by boundary
    cl_gdf = cl_gdf.to_crs(bd_gdf.crs)
    bd_geom = bd_gdf.loc[0, 'geometry']

    cl_in_bd = cl_gdf['geometry'].apply(bd_geom.intersection)
    if (rootDir):
        cl_in_bd.to_file(rootDir+'/cl_clip.shp')
    cl_in_bd.apply(cl_in_bd.intersection)

    ## start the algorithm from the longest segment
    idx = cl_in_bd.length.argmax()
    visited = {idx}
    stack = [cl_in_bd[idx]]

    def mergeAdj(seg, id):
        '''
        merge the segment with its adjacent
        '''
        visited.add(id)
        union_seg = seg.union(cl_in_bd[id])
        merged = ops.linemerge(union_seg)
        return merged

    ## find the longest merged lines
    candidate = []
    while (len(stack) > 0):
        seg = stack.pop()
        adj_ids = []
        for i in [0, -1]:
            end = Point(seg.coords[i])
            intc = cl_in_bd.intersection(end)
            intc = intc[~intc.is_empty]
            adj_id = [id for id in intc.index if id not in visited]
            if (len(adj_id) != 0):
                adj_ids.append(adj_id)
        if (len(adj_ids) == 0):
            candidate.append(seg)
        else:
            for id1 in adj_ids[0]:
                seg1 = mergeAdj(seg, id1)
                if (len(adj_ids) == 2):
                    for id2 in adj_ids[2]:
                        seg2 = mergeAdj(seg1, id2)
                        stack.append(seg2)
                else:
                    stack.append(seg1)
                
    ## select the one whose length is closet to boudary's approximated length
    approxL = get_approx_channel_length(bd_geom)
    candidate = gpd.GeoSeries(candidate)
    diff = np.abs(candidate.length - approxL)
    best_cl = candidate[np.argmin(diff)]
    cl2D = ops.transform(lambda x,y,z: tuple(filter(None, [x, y])), best_cl)
    cl_gdf = gpd.GeoDataFrame({'geometry': [cl2D]}, index=[0], crs=bd_gdf.crs)
    return cl_gdf
    

