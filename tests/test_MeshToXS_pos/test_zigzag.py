import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach

root_dir = os.path.dirname(__file__)
input_dir = f'{root_dir}/input/zigzag'
output_dir = f'{root_dir}/output/zigzag'

rr = RiverReach(
    centerline_file = f'{input_dir}/centerline.shp',
    boundary_file = f'{input_dir}/BoundaryPolygon.shp',
    mesh_file = f'{input_dir}/mesh_points.shp',
    survey_type = 'zigzag'
)

cst = CrossSectionTools(rr)

positions = [250, 300, 350, 400]
xs_file = f'{output_dir}/intp_xs.shp'
cst.MeshToXS(positions=positions, output_dir=output_dir, output_file=xs_file)