import os
from cstools.tools import CrossSectionTools
from cstools.reach import RiverReach

root_dir = os.path.dirname(__file__)
input_dir = f'{root_dir}/input/zigzag'
output_dir = f'{root_dir}/output/zigzag'

rr = RiverReach(
    centerline_file = f'{input_dir}/centerline.shp',
    boundary_file = f'{input_dir}/BoundaryPolygon.shp',
    observation_file = f'{input_dir}/observation.shp',
    survey_type = 'zigzag'
)

cst = CrossSectionTools(rr)

XS_file = f'{output_dir}/projected_XS.shp'
masks_file = f'{output_dir}/mask.shp'
mesh_point_file = f'{output_dir}/mesh_points.shp'
mesh_line_file = f'{output_dir}/mesh_lines.shp'

allpivots = cst.PointToMesh(with_original=True, 
                            output_dir=output_dir, 
                            output_point_file=mesh_point_file,
                            output_line_file=mesh_line_file,
                            num_vertices=21,
                            XS_file=XS_file,
                            masks_file=masks_file)