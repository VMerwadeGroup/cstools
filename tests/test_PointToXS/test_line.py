import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach


root_dir = os.path.dirname(__file__)
input_dir = f'{root_dir}/input/line'
output_dir = f'{root_dir}/output/line'

rr = RiverReach(
    centerline_file = f'{input_dir}/centerline.shp',
    boundary_file = f'{input_dir}/boundary.shp',
    observation_file = f'{input_dir}/observation.shp',
    survey_type = 'line'
)

cst = CrossSectionTools(rr)

XS_file = f'{output_dir}/projected_XS.shp'
masks_file = f'{output_dir}/mask.shp'
cst.PointToXS(output_dir=output_dir, XS_file=XS_file, masks_file=masks_file)