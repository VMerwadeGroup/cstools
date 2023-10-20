from __future__ import annotations

'''
    TODO List:
    1. Add some functions to check if the user provide enough data and information to launch the tools (especially the RasterTools)
    2. improve raster tools
    3. improve the calculation by JAX
'''
import os
import numpy as np
# import jax.numpy as jnp
from scipy import signal, interpolate
import pandas as pd
import geopandas as gpd
from shapely.geometry import *
from shapely import speedups
import warnings
import pyproj
from . import preprocess
from . import file_mangement

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
        output_dir = './temp'
        HD_ob = False
        shift = True
        smooth_cl = True
        unit_conversion_factor = 1 # length of unit : meter
        survey_type = 'zigzag' # default is zigzag (irregular shapes)
        cl_spacing = 10
        
        if 'centerline_file' in kwargs               :   centerline_file = kwargs.pop('centerline_file')
        if 'cross_section_file' in kwargs            :   cross_section_file = kwargs.pop('cross_section_file')
        if 'observation_file' in kwargs              :   observation_file = kwargs.pop('observation_file')
        if 'boundary_file' in kwargs                 :   boundary_file = kwargs.pop('boundary_file')
        if 'mesh_file' in kwargs                     :   mesh_file = kwargs.pop('mesh_file')
        if 'survey_type' in kwargs                   :   survey_type = kwargs.pop('survey_type')
        if 'map_crsid' in kwargs                     :   map_crsid = kwargs.pop('map_crsid')
        if 'output_dir' in kwargs                    :   output_dir = kwargs.pop('output_dir')
        if 'HD_ob' in kwargs                         :   HD_ob = kwargs.pop('HD_ob')
        if 'elev_base' in kwargs                     :   elev_base = kwargs.pop('elev_base')
        if 'shift' in kwargs                         :   shift = kwargs.pop('shift')
        if 'smooth_cl' in kwargs                     :   smooth_cl = kwargs.pop('smooth_cl')
        if 'unit_conversion_factor' in kwargs        :   unit_conversion_factor = kwargs.pop('unit_conversion_factor')
        if 'cl_spacing' in kwargs                    :   cl_spacing = kwargs.pop('cl_spacing')
      
        if (len(kwargs) != 0): raise ValueError('Please provide recognizable keyword arguments!')
        if (survey_type not in ['line', 'zigzag']): raise ValueError('Please assign the survey_type as "line" or "zigzag"!')
        
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

        boundary_driver = None
        if (observation_file):
            self.ob_gdf, self.lowest_elev = file_mangement.init_observation(
                                                path = observation_file,
                                                output_dir = output_dir,
                                                driver = None,
                                                use_HD = HD_ob,
                                                elev_base = elev_base,
                                                shift = shift
                                            )
            
            if (observation_file.startswith("https:") or observation_file.startswith("http:")):
                dir_name = os.path.basename(observation_file).replace(".ZIP", "")
                boundary_file = f"./{output_dir}/{dir_name}/{dir_name}.gdb"
                boundary_driver = "FileGDB"
                        
        if (boundary_file): 
            self.bd_gdf = file_mangement.init_boundary(
                    path = boundary_file,
                    output_dir = output_dir,
                    driver = boundary_driver
            )

        self.bd_geom = self.bd_gdf.geometry.iloc[0] if isinstance(self.bd_gdf, gpd.GeoDataFrame) else None
        self.user_crs = pyproj.CRS(map_crsid) if map_crsid else self.bd_gdf.crs
        self.proj_crs = self.user_crs
        self.checkCRS()

        if (cross_section_file):
            self.cs_gdf = file_mangement.init_cross_sections(path=cross_section_file)
        if (mesh_file):
            self.ms_gdf = file_mangement.init_mesh(path=mesh_file)
        if (centerline_file):
            self.cl_gdf = file_mangement.init_centerline(path=centerline_file, bd_file=self.bd_gdf)
        self.checkCRS()
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
            self.vertices_in_xs = int(self.approx_width / self.bw_org)
        if isinstance(self.cl_gdf, gpd.GeoDataFrame):
            ## ---------- centerline proporties for SN coordinates transform ----------
            self.coord_converter = preprocess.CoordConverter(self.cl_geom)
            self.cl_xyz = preprocess.get_xy_array_from_line(self.cl_geom)

    def checkCRS(self):
        if (self.proj_crs is None): raise ValueError(f'''Please assign a projected coordinate system or boudary layer!''')
        self.unit = self.proj_crs.coordinate_system.axis_list[1].unit_name

        ## check the target area is in which NAD 1983 UTM Zone
        if ((self.unit == 'degree')):
            centroid = self.bd_geom.centroid
            crs_zones = [pyproj.CRS.from_user_input(26910+i) for i in range(10)]
            for crs_zone in crs_zones:
                extent = crs_zone.area_of_use.bounds
                if (extent[0] < centroid.x < extent[2]) and (extent[1] < centroid.y < extent[3]):
                    self.proj_crs = crs_zone
                    break
            
        if (self.cl_gdf is not None):
            self.cl_gdf.to_crs(self.proj_crs, inplace=True)
            self.cl_geom = self.cl_gdf.loc[0, 'geometry']
        if (self.ob_gdf is not None): self.ob_gdf.to_crs(self.proj_crs, inplace=True)
        if (self.cs_gdf is not None): self.cs_gdf.to_crs(self.proj_crs, inplace=True)
        if (self.ms_gdf is not None): self.ms_gdf.to_crs(self.proj_crs, inplace=True)
        if (self.bd_gdf is not None): 
            self.bd_gdf.to_crs(self.proj_crs, inplace=True)
            self.bd_geom = self.bd_gdf.loc[0, 'geometry']

    def convert_sn_coord_for_layer(self, input = None, output_shp = None, output_csv = None, fields = None):
        '''
        Transfrom layer to S-N coordinates
        Input:
            input: path of file that can be read by geopandas or GeoDataFrame
            output_file: the path to store the shapefile
            output_csv: the path to export attributes as .csv or not
            fields: the fields will be copied from the original layer
        Output:
            save the shape file with S-N in attributes table
            save the .csv file (if output_csv exists)
        '''
        #print('Transforming to S-N coordinate...')
        
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

        # xyz_array = preprocess.get_xyz_array(gdf['geometry'])
        # sn_array = preprocess.convert_xy2sn(xyz_array[:, :2], self.cl_xyz, self.bw_org)
        
        # gdf['x'] = xyz_array[:, 0]
        # gdf['y'] = xyz_array[:, 1]
        # gdf['z'] = xyz_array[:, 2]
        # gdf['S'] = sn_array[:, 0]
        # gdf['N'] = sn_array[:, 1]

        df = pd.DataFrame(gdf.drop(columns=['geometry']))
        if output_shp:
            gdf.to_file(output_shp)
        
        if output_csv:
            df.to_csv(output_csv)

        return gdf
