from __future__ import annotations

import os
import urllib
import zipfile
import requests
from glob import glob
import geopandas as gpd
import fiona
from shapely.geometry import LineString
from ..preprocess.basic_features import *


def load_file_from_ehydro_url(url: str, root_dir : str ='./temp') -> str:
    survey_code = os.path.basename(url).replace(".ZIP", "")
    workdir = f'{root_dir}/{survey_code}'
    if (not os.path.isdir(workdir)):
        os.makedirs(workdir)
    else:
        ## TODO: have to make sure there are required files in workdir
        ## the file has been downloaded before
        return workdir
    
    ## the file did not be downloaded before
    ## download it from eHydro by the url
    try:   
        fname = os.path.basename(url)
        if (".ZIP" in fname):
            fname = fname.replace(".ZIP", ".zip")
        fname, _ = urllib.request.urlretrieve(url, filename=f'{root_dir}/{fname}')
        with zipfile.ZipFile(fname, 'r') as zip_ref:
            zip_ref.extractall(workdir)
        return workdir
    except:
        raise IOError(f'''Path "{url}" does not exist!, please give a valid file path or url''')

def load_flowline_from_rimorphis_api(url: str) -> gpd.GeoDataFrame:
    try:
        r = requests.get(url, verify=False) 
        j = r.json()
        coords=[]
        crs = 'EPSG:4326' ## assume the centerline from rimorphis is WGS 84 
        ## get the contents from json list
        if (isinstance(j, list)):
            j = j[0]
        if (isinstance(j, dict)):
            if ('features' in j):
                if (len(j['features']) > 1):
                    gdf = gpd.GeoDataFrame.from_features(j, crs=crs)
                    return gdf
                else:
                    coords = j['features'][0]['geometry']['coordinates'][0]
            else:
                coords = j['geometry']['coordinates'][0]                    
        else:
            raise IOError(f'''Path "{url}" does not return the required information! Object: {j}''')
        xy = [p[:2] for p in coords]
        cl_geom = LineString(xy)
        dic = {'id':[0], 'geometry': [cl_geom]}
        gdf = gpd.GeoDataFrame(dic, crs=crs)
        return gdf
    except:
        raise ValueError(f'''URL: "{url}" does not exist or not include required information!, please give a valid url.''')