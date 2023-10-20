from __future__ import annotations


import pandas as pd
import geopandas as gpd
import numpy as np
import os

class TrendModifier(object):
    def __init__(self, ob_file: str|pd.DataFrame|gpd.GeoDataFrame) -> None:
        df = pd.DataFrame()
        if (isinstance(ob_file, str)):
            if (os.path.exists(ob_file)):
                df = gpd.read_file(ob_file)
        elif (isinstance(ob_file, gpd.GeoDataFrame) or isinstance(ob_file, pd.DataFrame)):
            df = ob_file
        else:
            raise IOError("Please provide a valid file path, DataFrame, or GeoDataFrame.")
        self.ob_df = df
        self.poly_fit_s = None
        self.poly_fit_ng = None
        self.poly_fit_nl = None

    def remove_trend_in_s(self) -> None:
        p = np.polynomial.Polynomial.fit(self.ob_df['S'], self.ob_df['z'], deg=3, domain=None, window=None)
        z_trend = np.polynomial.polynomial.polyval(self.ob_df['S'], p.convert().coef)

        self.poly_fit_s = p
        self.ob_df['trend_s'] = z_trend
        self.ob_df['z_detrend_s'] = self.ob_df['z'] - z_trend

    def add_trend_in_s(self, ob_file: str|pd.DataFrame|gpd.GeoDataFrame) -> np.array:
        if (self.poly_fit_s is None):
            raise ValueError("Please fit the trend by orignal observations first.")
        
        df = pd.DataFrame()
        if (isinstance(ob_file, str)):
            if (os.path.exists(ob_file)):
                df = gpd.read_file(ob_file)
        elif (isinstance(ob_file, gpd.GeoDataFrame) or isinstance(ob_file, pd.DataFrame)):
            df = ob_file
        else:
            raise IOError("Please provide a valid file path, DataFrame, or GeoDataFrame.")
        df['s_trend'] = np.polynomial.polynomial.polyval(df['S'], self.poly_fit_s.convert().coef)

    def remove_trend_in_n_global(self) -> None:
        p = np.polynomial.Polynomial.fit(self.ob_df['N'], self.ob_df['z_detrend_s'], deg=3, domain=None, window=None)
        z_trend = np.polynomial.polynomial.polyval(self.ob_df['N'], p.convert().coef)

        self.poly_fit_ng = p
        self.ob_df['trend_ng'] = z_trend
        self.ob_df['z_detrend_ng'] = self.ob_df['z'] - z_trend
    
    def add_trend_in_n_global(self, ob_file: str|pd.DataFrame|gpd.GeoDataFrame) -> np.array:
        if (self.poly_fit_ng is None):
            raise ValueError("Please fit the trend by orignal observations first.")
        
        df = pd.DataFrame()
        if (isinstance(ob_file, str)):
            if (os.path.exists(ob_file)):
                df = gpd.read_file(ob_file)
        elif (isinstance(ob_file, gpd.GeoDataFrame) or isinstance(ob_file, pd.DataFrame)):
            df = ob_file
        else:
            raise IOError("Please provide a valid file path, DataFrame, or GeoDataFrame.")
        
        df['ng_trend'] = np.polynomial.polynomial.polyval(df['N'], self.poly_fit_ng.convert().coef)

    def remove_trend_in_n_local(self, window: float) -> None:
        
        half_w = window / 2
        smin, smax = self.ob_df['S'].round().min(), self.ob_df['S'].round().max()
        s_ls = np.arange(smin + half_w, smax + half_w, window)

        self.window = window
        self.s_ls = s_ls
        
        self.ob_df['z_detrend_nl'] = self.ob_df['z_detrend_s']
        self.ob_df['trend_nl'] = self.ob_df['z_detrend_s']
        self.poly_fit_nl = []
        # for s in range(int(smin), int(smax)+1):
        for s in s_ls:
            filt = (self.ob_df['S'] <= s + half_w) & (self.ob_df['S'] >= s - half_w)
            local_df = self.ob_df[filt]
            p = np.polynomial.Polynomial.fit(local_df['N'], local_df['z_detrend_s'], deg=3, domain=None, window=None)
            self.poly_fit_nl.append(p)
        
            z_trend = np.polynomial.polynomial.polyval(self.ob_df.loc[filt, 'N'], p.convert().coef)
            self.ob_df.loc[filt, 'trend_nl'] = z_trend
            self.ob_df.loc[filt, 'z_detrend_nl'] -= z_trend

    def add_trend_in_n_local(self, ob_file: str|pd.DataFrame|gpd.GeoDataFrame) -> np.array:
        if (self.poly_fit_nl is None):
            raise ValueError("Please fit the trend by orignal observations first.")
        
        df = pd.DataFrame()
        if (isinstance(ob_file, str)):
            if (os.path.exists(ob_file)):
                df = gpd.read_file(ob_file)
        elif (isinstance(ob_file, gpd.GeoDataFrame) or isinstance(ob_file, pd.DataFrame)):
            df = ob_file
        else:
            raise IOError("Please provide a valid file path, DataFrame, or GeoDataFrame.")
        
        half_w = self.window / 2
        for s, p in zip(self.s_ls, self.poly_fit_nl):
            filt = (df['S'] <= s + half_w) & (df['S'] >= s - half_w)
            local_df = df[filt]
            z_trend = np.polynomial.polynomial.polyval(local_df['N'], p.convert().coef)
            df.loc[filt, 'nl_trend'] = z_trend
    