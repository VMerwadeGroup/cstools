import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach

root_dir = os.path.dirname(__file__)
input_dir = f'{root_dir}/input/line'
output_dir = f'{root_dir}/output/line'

rr = RiverReach(
    centerline_file = f'{input_dir}/centerline.shp',
    boundary_file = f'{input_dir}/boundary.shp',
    mesh_file = f'{input_dir}/mesh_points.shp',
    survey_type = 'line',
    smooth_cl = False
)

cst = CrossSectionTools(rr)

positions = [400, 500, 750, 900]
xs_file = f'{output_dir}/intp_xs.shp'
cst.MeshToXS(positions=positions, output_dir=output_dir, output_file=xs_file)