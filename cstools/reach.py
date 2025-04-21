from __future__ import annotations

import os
# import jax.numpy as jnp
import pandas as pd
import geopandas as gpd
from shapely import speedups
import warnings
import pyproj
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Polygon, Point
from scipy import signal
from scipy import interpolate
import os
from . import preprocess
from . import file_mangement
from requests.exceptions import ConnectionError


gpd.options.use_pygeos = False

'''
PyGEOS is a C/Python library with vectorized geometry functions.
It could largely accelerate the computation of GeoDataFrame.
However, The use of PyGEOS is still experimental for geopandas.
Many unexpected errors happened, so I closed this library until it became stable.
'''

# import matplotlib.pyplot as plt
# from tqdm import tqdm

warnings.filterwarnings('ignore')
speedups.enable()

class RiverReach(object):
    def __init__(self, *args, **kwargs):
        if args:
            raise ValueError(  
                   'The cross-section tools constructor can only be called with keyword arguments for \
                    the following keywords: centerline_file, cross_section_file, boundary_file, observation_file, \
                    mesh_file, survey_type, map_crsid')
        
        ## default values for input parameters
        centerline_file = cross_section_file = boundary_file = observation_file = \
        mesh_file = survey_type = map_crsid = elev_base = None
        download_dir = './temp'
        HD_ob = False
        shift = True
        smooth_cl = True
        clip_by_boundary = False
        unit_conversion_factor = 1 # length of unit : meter
        survey_type = 'zigzag' # default is zigzag (irregular shapes)
        cl_spacing = 1
        input_driver = None
        positive_depth = True
        
        if 'centerline_file' in kwargs               :   centerline_file = kwargs.pop('centerline_file')
        if 'cross_section_file' in kwargs            :   cross_section_file = kwargs.pop('cross_section_file')
        if 'observation_file' in kwargs              :   observation_file = kwargs.pop('observation_file')
        if 'boundary_file' in kwargs                 :   boundary_file = kwargs.pop('boundary_file')
        if 'mesh_file' in kwargs                     :   mesh_file = kwargs.pop('mesh_file')
        if 'survey_type' in kwargs                   :   survey_type = kwargs.pop('survey_type')
        if 'map_crsid' in kwargs                     :   map_crsid = kwargs.pop('map_crsid')
        if 'download_dir' in kwargs                  :   download_dir = kwargs.pop('download_dir')
        if 'HD_ob' in kwargs                         :   HD_ob = kwargs.pop('HD_ob')
        if 'elev_base' in kwargs                     :   elev_base = kwargs.pop('elev_base')
        if 'shift' in kwargs                         :   shift = kwargs.pop('shift')
        if 'smooth_cl' in kwargs                     :   smooth_cl = kwargs.pop('smooth_cl')
        if 'unit_conversion_factor' in kwargs        :   unit_conversion_factor = kwargs.pop('unit_conversion_factor')
        if 'cl_spacing' in kwargs                    :   cl_spacing = kwargs.pop('cl_spacing')
        if 'input_driver' in kwargs                  :   input_driver = kwargs.pop('input_driver')
        if 'clip_by_boundary' in kwargs              :   clip_by_boundary = kwargs.pop('clip_by_boundary')
        if 'positive_depth' in kwargs                :   positive_depth = kwargs.pop('positive_depth')
      
        if (len(kwargs) != 0): raise ValueError(f'{kwargs} are not valid keyword arguments!')
        if (survey_type not in ['linear', 'zigzag']): raise ValueError('Please assign the survey_type as "linear" or "zigzag"!')
        
        ## update settings based on the input parameters
        self.HD_ob = HD_ob
        self.elev_base = elev_base
        self.shift = shift
        self.lowest_elev = None
        self.unit = None
        self.smooth_cl = smooth_cl
        self.unit_conversion_factor = unit_conversion_factor

        self.survey_type = survey_type
        self.cl_gdf = self.cs_gdf = self.ob_gdf = self.bd_gdf = self.ms_gdf = None

        ## getting files from different methods
        ## 1. eHydro_url
        ## 2. eHydro_folder
        ## 3. .gdb, TODO: get the layer name from input   
        ## 4. file path

        if (observation_file):
            self.ob_gdf, self.lowest_elev = file_mangement.init_observation(
                                                path = observation_file,
                                                download_dir = download_dir,
                                                driver = input_driver,
                                                use_HD = HD_ob,
                                                elev_base = elev_base,
                                                shift = shift,
                                                positive_depth=positive_depth
                                            )
            
            if (input_driver is not None):
                if ('eHydro' in input_driver) and (boundary_file is None):
                    driver = "FileGDB"
                    layer = "SurveyJob"
                    if (input_driver == 'eHydro_url'):
                        site_ID = os.path.basename(observation_file).replace(".ZIP", "")
                        boundary_file = os.path.join(download_dir, site_ID, f"{site_ID}.gdb")
                    elif (input_driver == 'eHydro_dir'):
                        site_ID = os.path.basename(observation_file)
                        boundary_file = os.path.join(observation_file, f"{site_ID}.gdb")
                    else:
                        raise ValueError("Invalid value of 'input_driver'.")
                    self.bd_gdf = file_mangement.init_boundary(path = boundary_file, driver = driver, layer = layer)

        ## TODO: check if boundary file should from .gdb or not        
        if (boundary_file) and (self.bd_gdf is None): 
            self.bd_gdf = file_mangement.init_boundary(path = boundary_file)

        self.bd_geom = self.bd_gdf.geometry.iloc[0] if isinstance(self.bd_gdf, gpd.GeoDataFrame) else None
        self.user_crs = pyproj.CRS(map_crsid) if map_crsid else self.bd_gdf.crs
        self.proj_crs = self.user_crs
        self._checkCRS()

        ## only process data in the boundary
        if (clip_by_boundary):
            if (self.ob_gdf is not None) and (self.bd_geom is not None):
                self.ob_gdf = self.ob_gdf[self.ob_gdf.within(self.bd_geom)]
                self.ob_gdf = self.ob_gdf.reset_index(drop=True)
            else:
                raise ValueError("The observation or boundary layers are not initialized correctly!")

        if (cross_section_file):
            self.cs_gdf = file_mangement.init_cross_sections(path=cross_section_file)
        if (mesh_file):
            self.ms_gdf = file_mangement.init_mesh(path=mesh_file)
        if (centerline_file):
            self.cl_gdf = file_mangement.init_centerline(path=centerline_file, bd_file=self.bd_gdf)
        else:
            cl_geom = preprocess.create_cl_from_bd(self.bd_geom)
            self.cl_gdf = gpd.GeoDataFrame(geometry=[cl_geom], index=[0], crs=self.proj_crs)
        self._checkCRS()
        ## make sure the all CRS are the same and projected; then start processing geometries    
       
        ## ---------- defalut values for parameters ----------
        ## unit_converion_factor: length of unit / meter
        ## bw_org is 5 meters, the value will adjust with the unit
        self.unit_conversion_factor = self.proj_crs.coordinate_system.axis_list[1].unit_conversion_factor 
        self.bw_org = 5 / self.unit_conversion_factor
        if (smooth_cl):        
            self.cl_gdf.loc[0, 'geometry'] = preprocess.smooth_line_geom(self.cl_gdf.loc[0, 'geometry'], spacing = cl_spacing)
        self.cl_geom = self.cl_gdf.loc[0, 'geometry'] if isinstance(self.cl_gdf, gpd.GeoDataFrame) else None

        if (isinstance(self.cl_gdf, gpd.GeoDataFrame) and isinstance(self.bd_gdf, gpd.GeoDataFrame)):
            self.approx_width = preprocess.get_approx_channel_width(bd_geom=self.bd_geom)
            self.interp_xs_space = self.approx_width / 5
            self.vertices_in_xs = max(int(self.approx_width / self.bw_org), 21)
        if isinstance(self.cl_gdf, gpd.GeoDataFrame):
            ## ---------- centerline proporties for SN coordinates transform ----------
            self.coord_converter = preprocess.CoordConverter(self.cl_geom)
            self.cl_xyz = preprocess.get_xy_array_from_line(self.cl_geom)

    def _checkCRS(self):
        if (self.proj_crs is None): 
            raise ValueError("Please assign a projected coordinate system or a boudary layer!")
        self.unit = self.proj_crs.coordinate_system.axis_list[1].unit_name

        ## check the target area is in which NAD 1983 UTM Zone in CONUS, if no projected CRS is provided
        ## The default CRS is NAD 1983 UTM Zone 10N to 19N in meters
        if (self.unit == 'degree'):
            self.proj_crs = self.bd_gdf.estimate_utm_crs()
            self.unit = self.proj_crs.coordinate_system.axis_list[1].unit_name
            
        ## update the CRS in every layer and geometry to make them consistent
        if (self.ob_gdf is not None): self.ob_gdf.to_crs(self.proj_crs, inplace=True)
        if (self.cs_gdf is not None): self.cs_gdf.to_crs(self.proj_crs, inplace=True)
        if (self.ms_gdf is not None): self.ms_gdf.to_crs(self.proj_crs, inplace=True)
        if (self.cl_gdf is not None):
            self.cl_gdf.to_crs(self.proj_crs, inplace=True)
            self.cl_geom = self.cl_gdf.loc[0, 'geometry']
        if (self.bd_gdf is not None): 
            self.bd_gdf.to_crs(self.proj_crs, inplace=True)
            self.bd_geom = self.bd_gdf.loc[0, 'geometry']

    def convert_sn_coord_for_layer(self, 
                                   input: str|gpd.GeoDataFrame|None = None, 
                                   output_shp: str|None = None, 
                                   output_csv: str|None = None
                                   ) -> gpd.GeoDataFrame:
        """Convert the point coordinates from Cartiesian system to the curvilier system.
        The coordinates are converted based on the centerline of this object.

        Parameters
        ----------
        input : str | gpd.GeoDataFrame | None, optional
            The input file of the observation points. 
            By default None, the ob_gdf in this object will be used.
            It can also be an existed GeoDataFrame, or a path of an importable file for geopandas.
        output_shp : str | None, optional
            The path of a output georeferenced file, by default None
        output_csv : str | None, optional
            The path of a output csv file, by default None

        Returns
        -------
        gpd.GeoDataFrame
            A GoeDataFrame including x, y, z, S, N in its attribute table.

        Raises
        ------
        IOError
            Input in invalid.
        """
        
        if (input is None):
            gdf = self.ob_gdf
        elif (isinstance(input, gpd.GeoDataFrame)):
            gdf = input
            gdf.to_crs(self.proj_crs)
        elif (os.path.exists(input)):
            gdf = gpd.read_file(input)
            gdf = gdf.to_crs(self.proj_crs)
        else:
            raise IOError("Please provide a GeoDataFrame or valid file path!")
            
        getXYZFromGDF = lambda row: (row.geometry.coords[0])
        XYZ = gdf.apply(getXYZFromGDF, axis=1)

        SN = XYZ.apply(lambda row: self.coord_converter.xy2sn_coord(*row))
        (x, y, z) = zip(*XYZ.values)
        s, n, _ = zip(*SN)

        gdf['x'] = x
        gdf['y'] = y
        gdf['z'] = z
        gdf['S'] = s
        gdf['N'] = n

        df = pd.DataFrame(gdf.drop(columns=['geometry']))
        if output_shp:
            gdf.to_file(output_shp)
        
        if output_csv:
            df.to_csv(output_csv)

        return gdf
