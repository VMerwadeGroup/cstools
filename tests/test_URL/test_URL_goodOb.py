import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach


root_dir = os.path.dirname(__file__)
output_dir = f'{root_dir}/output/goodOb'

rr = RiverReach(
    centerline_file = 'https://rimorphis.org/platform/api/postgrest/river_flowline?and=(foreign_id.eq.22000600016702)&select=geometry',
    observation_file = 'https://ehydrotest.blob.core.usgovcloudapi.net/ehydro-surveys/CEMVP/UM_SP_P02_20190723_CS_8354_8364.ZIP',
    survey_type = 'zigzag'
)

cst = CrossSectionTools(rr)

XS_file = f'{output_dir}/projected_XS.shp'
mesh_point_file = f'{output_dir}/mesh_points.shp'
mesh_line_file = f'{output_dir}/mesh_lines.shp'

cst.PointToXS(output_dir=output_dir, XS_file=XS_file)
allpivots = cst.PointToMesh(with_original=True, output_dir=output_dir, 
                            output_point_file=mesh_point_file,
                            output_line_file=mesh_line_file)

# %%
