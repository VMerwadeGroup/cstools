import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach

root_dir = os.path.dirname(__file__)
input_dir = f'{root_dir}/input/zigzag'
output_dir = f'{root_dir}/output/zigzag'

rr = RiverReach(
    centerline_file = f'{input_dir}/centerline.shp',
    boundary_file = f'{input_dir}/BoundaryPolygon.shp',
    cross_section_file = f'{input_dir}/projected_XS.shp',
    survey_type = 'zigzag'
)

cst = CrossSectionTools(rr)

mesh_point_file = f'{output_dir}/mesh_points.shp'
mesh_line_file = f'{output_dir}/mesh_lines.shp'
allpivots = cst.XSToMesh(with_original=True, output_dir=output_dir, 
                         output_point_file=mesh_point_file,
                         output_line_file=mesh_line_file)