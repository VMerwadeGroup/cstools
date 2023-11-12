import os
from cstools.reach import RiverReach
import geopandas as gpd
from shapely.geometry import Polygon

os.chdir(os.path.dirname(__file__))
import pytest

def test_eHydro_dir():
    rr = RiverReach(
        observation_file = "./LM_02_GCG_20180906_CS_3454_3660",
        centerline_file = "centerline.geojson",
        survey_type = 'zigzag',
        input_driver = 'eHydro_dir'
    )
    
    ## check observation
    if (not isinstance(rr.ob_gdf, gpd.GeoDataFrame)):
        assert False, "The observation GeoDataFrame does not exist!"    
    if (len(rr.ob_gdf) < 10):
        assert False, "The observation GeoDataFrame does not load correct data!"

    ## check boundary
    if (not isinstance(rr.bd_gdf, gpd.GeoDataFrame)):
        assert False, "The boundary GeoDataFrame does not exist!"    
    if (not isinstance(rr.bd_geom,  Polygon)):
        assert False, "The boundary geometry is not a Polygon!"

    ## check CRS
    if (rr.ob_gdf.crs != rr.bd_gdf.crs):
        assert False, "The coordinate reference systems (CRS) should be consistent between observation and boundary!"

    assert True
