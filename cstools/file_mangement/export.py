from __future__ import annotations

import os
import geopandas as gpd

def check_dir_exists(dir_path: str) -> None:
    if not os.path.exists(dir_path):
        raise IOError(f'''Directory "{dir_path}" does not exist!''')

def get_output_path(dir: str, fname: str|None) -> str|None:
    if (fname is None):
        return None
    elif dir:
        return os.path.join(dir, os.path.basename(fname))
    else:
        return fname

def export_geofile(gdf: gpd.GeoDataFrame, output_file: str) -> None:
    output_dir = os.path.realpath(os.path.dirname(output_file))
    if (not os.path.exists(output_dir)):
        os.makedirs(output_dir)
    
    gdf.to_file(output_file)
    
def __export_geofile_rimorphis(gdf: gpd.GeoDataFrame, output_file: str) -> None:
    output_root_dir = os.path.realpath(os.path.dirname(output_file))
    output_dirs = [f'{output_root_dir}/shp', f'{output_root_dir}/geojson', f'{output_root_dir}/geojson_4326']
    for d in output_dirs:
        if (not os.path.exists(d)):
            os.makedirs(d)

    file_name = os.path.basename(output_file)
    file_name_ls = file_name.split('.')
    name, ext = file_name_ls[0], file_name_ls[1]
    types = ['.shp', '.geojson', '.geojson_4326']
    if ext not in types:
        types.append(f'.{ext}')
    for d, t in zip(output_dirs, types):
        if t=='.geojson_4326':
            ## store geojson file in EPSG:4326 for visulization
            ## filename: name.geojson
            gdf_4326 = gdf.to_crs(4326)
            gdf_4326.to_file(f'{d}/{name}{".geojson"}', driver='GeoJSON')
        elif t=='.geojson':
            gdf.to_file(f'{d}/{name}{t}', driver='GeoJSON')
        else:
            gdf.to_file(f'{d}/{name}{t}')