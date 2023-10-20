import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach

root_dir = os.path.dirname(__file__)
output_dir = f'{root_dir}/output/multiCl'

rr = RiverReach(
	centerline_file = 'https://rimorphis.org/platform/api/postgrest/rpc/get_flowlines_in_bb?lngw=-93.0703916041091&lats=44.94795430845617&lnge=-93.07239160410911&latn=44.94595430845617',
	observation_file = 'https://ehydrotest.blob.core.usgovcloudapi.net/ehydro-surveys/CEMVP/UM_SP_P02_20211006_CS_8380_8388.ZIP',
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