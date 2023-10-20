from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from skgstat import DirectionalVariogram, Variogram
from skgstat import OrdinaryKriging
import matplotlib.pyplot as plt
import seaborn as sns
import os
from .preprocess import CoordConverter
from .reach import *
from joblib import Parallel, delayed
import multiprocessing as mp
from functools import partial

class AnisotropicKriging(object):
    def __init__(self) -> None:
        self.anisotropy_factor = 1.
        self.v = None
        self.ok = None
        self.sigma = None        

    def get_anisotropy_factor(self, ob_file: str|pd.DataFrame|gpd.GeoDataFrame) -> float:
        df = pd.DataFrame()
        if (isinstance(ob_file, str)):
            if (os.path.exists(ob_file)):
                df = gpd.read_file(ob_file)
        elif (isinstance(ob_file, gpd.GeoDataFrame) or isinstance(ob_file, pd.DataFrame)):
            df = ob_file
        else:
            raise IOError("Please provide a valid file path, DataFrame, or GeoDataFrame.")
        
        self.dv_s = DirectionalVariogram(df[['S', 'N']].values, 
                                         df['z'].values, 
                                         azimuth=0,
                                         tolerance=45,
                                         bandwidth=50,
                                         maxlag=0.8, 
                                         n_lags=25, 
                                         model='spherical', 
                                         directional_model='triangle',
                                         normalize=False)
        res_s = self.dv_s.describe()
        del self.dv_s

        self.dv_n = DirectionalVariogram(df[['S', 'N']].values,
                                         df['z'].values, 
                                         azimuth=90,
                                         tolerance=45,
                                         bandwidth=50,
                                         maxlag=0.8, 
                                         n_lags=25, 
                                         model='spherical', 
                                         directional_model='triangle',
                                         normalize=False)
        res_n = self.dv_n.describe()
        del self.dv_n

        aniso = res_n['effective_range'] / res_s['effective_range']
        self.anisotropy_factor = aniso
        return aniso

    def interp(self,
               ob_file: str|pd.DataFrame|gpd.GeoDataFrame,
               s_spacing: float = 2,
               n_spacing: float = 2
               ) -> tuple[np.array, np.array, np.array]:
        
        df = pd.DataFrame()
        if (isinstance(ob_file, str)):
            if (os.path.exists(ob_file)):
                df = gpd.read_file(ob_file)
        elif (isinstance(ob_file, gpd.GeoDataFrame) or isinstance(ob_file, pd.DataFrame)):
            df = ob_file
        else:
            raise IOError("Please provide a valid file path or GeoDataFrame.")

        z_val = df['z'].values

        anisotropy_r = self.get_anisotropy_factor(ob_file=ob_file)
        df['S_re'] = df['S'] * anisotropy_r

        self.v = Variogram(df[['S_re', 'N']].values, 
                           z_val, 
                           maxlag=0.8, 
                           n_lags=25, 
                           model='spherical')
        
        self.ok = OrdinaryKriging(self.v, min_points=5, max_points=15)

        s_spacing_re = s_spacing * anisotropy_r
        s_re = np.arange(df['S_re'].min(), df['S_re'].max()+s_spacing_re, s_spacing_re)
        n = np.arange(df['N'].min(), df['N'].max()+n_spacing, n_spacing)

        ss, nn = np.meshgrid(s_re, n)
        field = self.ok.transform(ss.flatten(), nn.flatten()).reshape(ss.shape)

        ### TODO: find the way to fix that kriging provides some outliers
        z_max = df['z'].max()
        z_min = df['z'].min()
        field[field > z_max] = z_max
        field[field < z_min] = z_min
        
        # self.sigma = self.ok.sigma.reshape(ss.shape)
        del self.v
        del self.ok

        s = s_re / anisotropy_r

        return s, n, field
    
    def get_grids_in_boundary(self, s:np.array, n:np.array, z:np.array,
                              reach:RiverReach) -> gpd.GeoDataFrame:
        
        ss, nn = np.meshgrid(s, n)
        df = pd.DataFrame({'S':ss.flatten(), 'N':nn.flatten(), 'z':z.flatten()})
        df = df.sort_values(by=['S', 'N'])
        df = df.reset_index(drop=True)

        df_sn2xy = lambda row: Point(*reach.coord_converter.sn2xy_coord(*row))
        df['geometry'] = df[['S', 'N', 'z']].apply(df_sn2xy, axis=1)

        # gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=reach.proj_crs)
        # gdf = gdf.clip(reach.bd_geom, keep_geom_type=True)
        # gdf = gdf.sort_values(by=['S', 'N'])
        # gdf = gdf.reset_index(drop=True)
        return df


