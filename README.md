# cstools

cstools is a Python package that provides functionality for creating curvilinear mesh through interpolations.


## Installation

You can install cstools using ``pip``:

    $ pip install cstools

or using ``conda`` (``mamba``):

    $ conda install -c conda-forge cstools


## Quick start

The following example shows how to utilize this pacakge to create curvilinear points and a denser mesh.

A ``RiverReach`` class in ``cstools.reach`` is used to store the information of a river reach. The required inputs include ``centerline_file``, ``boundary_file``, and ``observation_file``. A function, ``convert_sn_coord_for_layer``, in the RiverReach object can be used to convert coordinates from the Cartician system (x, y) to cuvilinear system (S, N).

A ``CrossSectionTools`` in ``cstools.tools`` is used to interpolate points or cross-sections for mesh constructions. The default parameters are given or will be calculated when the functions are called. Users also can customize their own parameters. The example here shows the argments to set up paths for output files.

```python

    from cstools.reach import RiverReach
    from cstools.tools import CrossSectionTools

    input_dir = '/input/zigzag'
    output_dir = '/output/zigzag'

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
```
