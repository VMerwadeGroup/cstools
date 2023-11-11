import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach

root_dir = os.path.dirname(__file__)
output_dir = f'{root_dir}/output/lackOb'

rr = RiverReach(
    centerline_file = 'https://rimorphis.org/platform/api/postgrest/river_flowline?and=(foreign_id.eq.22000600044323)&select=geometry',
    observation_file = 'https://ehydrotest.blob.core.usgovcloudapi.net/ehydro-surveys/CEMVP/UM_SP_SAF_20150316_CS_8554_8558.ZIP',
    survey_type = 'zigzag',
    input_driver = 'eHydro_url'
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
