from __future__ import annotations

import os
import geopandas as gpd
import pandas as pd
import fiona
from shapely.geometry import Point
from .load import *
from ..preprocess.basic_features import get_cl_from_multilines

'''
TODO:
1. initialize observation by tables
2. specified z by the column of the table, so that the 2D points will be allowed
'''


'''
if the file path is an exist file, directly access it
Now the URL is only supported observations and boudarys file from eHydro
The centerline from NHD+ (Iowa server)
'''

def init_centerline(path: str, bd_file: str|None = None) -> gpd.GeoDataFrame:
    if (path.startswith("https:") or path.startswith("http:")):
        gdf = load_flowline_from_rimorphis_api(path)
    elif (os.path.exists(path)):
        gdf = gpd.read_file(path)
    else:
        raise ValueError(f'The path/URL "{path}" does not exist')
    
    gdf = gdf.explode()
    gdf.reset_index(drop=True, inplace=True)
    if (len(gdf) > 1):
        if (not bd_file is None):
            if isinstance(bd_file, str):
                bd_gdf = gpd.GeoDataFrame(bd_file)
            elif isinstance(bd_file, gpd.GeoDataFrame):
                bd_gdf = bd_file
            gdf = get_cl_from_multilines(cl_gdf=gdf, bd_gdf=bd_gdf)
        else:
            raise ValueError("A boundary file is necessary to process multiple flowlines.")
    return gdf

def init_boundary(path: str, 
                  driver: str|None = None,
                  layer: str|None = None,
                  ) -> gpd.GeoDataFrame:
    if (driver is None):
        if (os.path.exists(path)):
            gdf = gpd.read_file(path)
        else:
            raise ValueError("Invalid path.")        
    elif (driver == "FileGDB"):
        ## default the .gdb file is from eHydro
        ## TODO: it should be updated for general input source
        gdf = gpd.read_file(path, driver=driver, layer=layer)
    else:
        raise ValueError("Invalid driver.")
    
    gdf = gdf.explode()
    gdf.reset_index(drop=True, inplace=True)

    if (len(gdf) > 1):
        raise ValueError("This tool only accept one boundary geometry.")
    return gdf

def init_mesh(path: str) -> gpd.GeoDataFrame:
    return gpd.read_file(path)

def init_cross_sections(path: str) -> gpd.GeoDataFrame:
    return gpd.read_file(path)

def init_observation(
        path: str,
        download_dir: str = './temp',
        driver: str|None = None,
        use_HD: bool = True,
        elev_base: float|bool|None = False,
        positive_depth: bool = True,
        shift: bool = True
    ) -> tuple[gpd.GeoDataFrame, float|None]:
    '''
    elev_base is True:
        If the observation data is water depth, the user should specify the elevation base to convert depth to elevations.
        The further computation is based on the elevation.
    shift is True:
        The shift is used for visualizaion.
        Since the elevations could be a large values, this flag can shift elevation values by subtracting the lowest elevation value.
        The origin of the elevation-axis would be set at the lowest elevation value in the whole dataset.
    '''

    elev_base = None if (elev_base is True) else elev_base
    lowest_elev = None

    if (driver is None):
        if (os.path.exists(path)):
            gdf = gpd.read_file(path)
            if (elev_base):
                lowest_elev = convert_depth_2_elev(gdf=gdf, elev_base=elev_base, positive_depth=positive_depth, shift=shift)
            elif (shift):
                lowest_elev = shift_elev_by_lowest(gdf=gdf)
        else:
            try:
                workdir = load_file_from_ehydro_url(url=path, root_dir=download_dir) 
                gdf, lowest_elev = process_xyz_in_workdir(workdir=workdir, use_HD=use_HD, elev_base=True, shift=shift)
            except:
                raise ValueError("Invalid path.")        
    elif (driver == "eHydro_url"):
        try:
            requests.get(path)
        except ConnectionError:
            print("Invalid url!")
        workdir = load_file_from_ehydro_url(url=path, root_dir=download_dir)
        gdf, lowest_elev = process_xyz_in_workdir(workdir=workdir, use_HD=use_HD, elev_base=True, shift=shift)
    elif (driver == 'eHydro_dir'):
        ## TODO: check if the directory has the necesary data
        if (not os.path.isdir(path)):
            raise ValueError("Invalid directory path of the eHydro dataset!")
        gdf, lowest_elev = process_xyz_in_workdir(workdir=path, use_HD=use_HD, elev_base=True, shift=shift)
    else:
        raise ValueError("Invalid driver.")
        
    gdf = gdf.explode()
    gdf.reset_index(drop=True, inplace=True)

    return gdf, lowest_elev


make_point_by_xyz = lambda row: Point(*row)

def find_elev_base_from_metadata_ehydro(workdir: str) -> float | None:
    elev_base = None
    XYZH_file = glob(f'{workdir}/*.XYZH')
    XYZ_file = glob(f'{workdir}/*.XYZ')

    ## try .XYZH file
    if (len(XYZH_file) > 0):
        with open(XYZH_file[0], 'r') as f:
            for line in f:
                if ('Water_Elevation' in line):
                    num = line.split('==')[1]
                    return float(num.replace('\n', ''))
    
    ## try .XYZ file
    elif (len(XYZ_file) > 0):
        with open(XYZ_file[0], 'r') as f:
            for line in f:
                if ('Water Surface Elevation' in line):
                    num = line.split(' ')[5]
                    elev_base = float(num.replace('ft', ''))
                if (line[0] != "#"):
                    break
    else:
        raise ValueError("The required files do not exist.")

    return elev_base

def process_xyz_in_workdir(
        workdir: str,
        use_HD: bool = True,
        elev_base: float|bool|None = False,
        positive_depth: bool = True,
        shift: bool = True
    ) -> tuple[gpd.GeoDataFrame, float|None]:

    '''
    positive depth: below water is positive, above is negative
    '''

    ## check if elev_base is provided
    if (elev_base is None) and (workdir):
        elev_base = find_elev_base_from_metadata_ehydro(workdir=workdir)
    
    if (elev_base is None):
        raise ValueError("Please provide valid elevation base.")

    survey_code = os.path.basename(workdir)
    crs = fiona.open(f"{workdir}/{survey_code}.gdb").crs
    xyz_A_file = f'{workdir}/{survey_code}_A.XYZ'
    xyz_file = f'{workdir}/{survey_code}.XYZ'

    ## use HD survey points, or preprocessed survey points
    if (use_HD):
        if (os.path.exists(xyz_A_file)):
            df = pd.read_csv(xyz_A_file, sep=" ", header=None, names=['x', 'y', 'z'])
        else:
            print('''Warning: Cannot find HD observations. Use preprocessed observation instead!''')
            df = pd.read_csv(xyz_file, sep=" ", header=None, names=['x', 'y', 'z'])
    else:
        df = pd.read_csv(xyz_A_file, sep=" ", header=None, names=['x', 'y', 'z'])

    ## if values represents depths, is downward a positive or negative direction
    if (elev_base):
        if (positive_depth):
            df['elev'] = elev_base - df['z']
        else:
            df['elev'] = elev_base + df['z']

    lowest_elev = None
    if (shift):
        lowest_elev = df['elev'].min()
        df['elev'] -= lowest_elev

    df['geometry'] = df[['x', 'y', 'elev']].apply(make_point_by_xyz, axis=1)
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=crs)

    return gdf, lowest_elev


def convert_depth_2_elev(
        gdf: gpd.GeoDataFrame, 
        elev_base: float,
        positive_depth: bool = False,
        shift: bool = True
    ) -> float|None:

    get_coords = lambda geom: geom.coords[0]
    coords = gdf['geometry'].apply(get_coords)
    ## there must be xyz in coords
    x, y, z = zip(*coords)
    gdf['x'] = x
    gdf['y'] = y
    gdf['z'] = z
    if (positive_depth):
        gdf['elev'] = elev_base - gdf['z']
    else:
        gdf['elev'] = elev_base + gdf['z']

    lowest_elev = None
    if (shift):
        lowest_elev = gdf['elev'].min()
        gdf['elev'] -= lowest_elev
    
    gdf['geometry'] = gdf[['x', 'y', 'elev']].apply(make_point_by_xyz, axis=1)

    return lowest_elev

def shift_elev_by_lowest(gdf: gpd.GeoDataFrame) -> float:

    get_coords = lambda geom: geom.coords[0]
    coords = gdf['geometry'].apply(get_coords)
    ## there must be xyz in coords
    x, y, z = zip(*coords)
    gdf['x'] = x
    gdf['y'] = y
    gdf['z'] = z

    lowest_elev = gdf['z'].min()
    gdf['z'] -= lowest_elev

    gdf['geometry'] = gdf[['x', 'y', 'z']].apply(make_point_by_xyz, axis=1)
    return lowest_elev




        