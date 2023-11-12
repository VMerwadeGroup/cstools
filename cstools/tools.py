from __future__ import annotations


from .reach import *
from .locate_cs import *
from . import curvilinear_interp
from . import file_mangement


class CrossSectionTools(object):
    """
    There four different methods for users.
    1. XSToMesh
    2. PointToXS
    3. PointToMesh
    4. MeshToXS

    The methods do not require any input from users.
    Users can use them without any settings.
    The inner functions will automatically determine the values of parameters according to the variables of the object.
    """

    def __init__(self, reach: RiverReach, output_dir: str='./temp'):
        self.reach = reach
        self.reach.output_dir = output_dir

    def update_cs_gdf(self, cs_df: pd.DataFrame, cs_file: str|None) -> None:
        getXSLineString = lambda xsgroup: LineString([Point(row[['x', 'y', 'z']]) for _, row in xsgroup.iterrows()])
        g = cs_df.groupby('XS_id')
        lines = g.apply(getXSLineString)
        self.reach.cs_gdf = gpd.GeoDataFrame({'geometry': lines}, crs=self.reach.cs_gdf.crs)
        if (cs_file):
            file_mangement.export_geofile(self.reach.cs_gdf, cs_file)

    def XSToMesh(self, num_vertices : int|None = None,
                 space_min : float|None = None,
                 space_max : float|None = None,
                 with_original: bool|None = False,
                 output_dir: str|None = None,
                 output_point_file: str|None = None,
                 output_line_file: str|None = None
                ) -> pd.DataFrame:
        
        """This function is used to generate mesh files (points or line strings) from cross-sections.

        Parameters
        ----------
        num_vertices : int, optional
            The number of vertices in each cross-section, by default None
        space_min : float, optional
            The min spacing between cross-sections, by default None
        space_max : float, optional
            The max spacing between cross-sections, by default None
        with_original : bool, optional
            To extract values from observed points after interpolating cross-sections or not, by default False
        output_dir : str, optional
            The path of directory to store the output files, by default None
        output_point_file : str, optional
            The path to save the mesh files as line features. by default None
        output_line_file : str, optional
            The path to save the mesh files as linestring features. by default None

        Returns
        -------
        DataFrame
            The DataFrame stores the x, y, S, N, and z coordinates of vertices in each cross-section.
            The XS_id in the table indicates which cross-section each point belongs to.
        """

        output_dir = output_dir if output_dir else self.reach.output_dir
        file_mangement.check_dir_exists(output_dir)

        output_point_file = file_mangement.get_output_path(output_dir, output_point_file)            
        output_line_file = file_mangement.get_output_path(output_dir, output_point_file)         

        self.reach.cl_gdf.to_crs(self.reach.proj_crs, inplace=True)
        self.reach.bd_gdf.to_crs(self.reach.proj_crs, inplace=True)
        self.reach.cs_gdf.to_crs(self.reach.proj_crs, inplace=True)
        with_original = False if self.reach.ob_gdf is None else with_original
        xs_df = get_cs_df_from_gdf(self.reach)
        allvertices = curvilinear_interp.generate_mesh(reach=self.reach, num_vertices=num_vertices, xs_df=xs_df, 
                                                       space_min=space_min, space_max=space_max, with_original=with_original, 
                                                       output_point_file=output_point_file,
                                                       output_line_file=output_line_file)
        return allvertices

    def PointToXS(self,
                  bins: int|None = None, 
                  num_vertices: int|None = None, 
                  positions: list[float]|None = None, 
                  bw: list[float]|float|None = None, 
                  dist_threshold: float = 2.0, 
                  width_threshold: float = 0.0, 
                  extend_to_boundary: bool = True,
                  output_dir: str|None = None,
                  XS_file: str|None = None, 
                  masks_file: str|None = None
                  ) -> pd.DataFrame:
        """This function is used to generate cross-section files (line strings) from a observed points file.

        Parameters
        ----------
        bins : int | None, optional
            The number of bins for the signal function to find out the locations of cross-sections, by default None
        num_vertices : int | None, optional
            The number of vertices in each cross-section, by default None
        positions : list[float] | None, optional
            Given the positions of cross-sections by S-coordinates, by default None
        bw : list[float] | float | None, optional
            The width of buffers that could be set as one value (float) or as different values (list), by default None
        dist_threshold : float, optional
            The maximum distance between two cross-sections.
            If two cross-sections are larger than this threshold, a cross-section will be inserted to the middel, by default 2.0
        width_threshold : float, optional
            The ratio of the cross-section length and channel width should larger than this threshold, or it will be removed, by default 0.0
        extend_to_boundary : bool, optional
            To extend the interpolation range to the boundary or not, by default True
        output_dir : str | None, optional
            The path of directory to store the output files, by default None
        XS_file : str | None, optional
            The file name or path of the cross-section file.
            If output_dir is specified, it will be saved in output_dir, by default None
        masks_file : str | None, optional
            The file name or path of the file that is used to clip points for cross-sections.
            If output_dir is specified, it will be saved in output_dir, by default None

        Returns
        -------
        pd.DataFrame
            The DataFrame stores the x, y, S, N, and z coordinates of vertices in each cross-section.
            The XS_id in the table indicates which cross-section each point belongs to.
        """
        
        output_dir = output_dir if output_dir else self.reach.output_dir
        file_mangement.check_dir_exists(output_dir)
        XS_file = file_mangement.get_output_path(output_dir, XS_file)
        masks_file = file_mangement.get_output_path(output_dir, masks_file)
        
        self.reach.ob_gdf.to_crs(self.reach.proj_crs, inplace=True)
        self.reach.cl_gdf.to_crs(self.reach.proj_crs, inplace=True)
        self.reach.bd_gdf.to_crs(self.reach.proj_crs, inplace=True)
        
        self.ob_gdf = self.reach.convert_sn_coord_for_layer()
        locate_cs_by_hist(self.reach, bins, num_vertices, positions, bw, extend_to_boundary, dist_threshold=dist_threshold, width_threshold=width_threshold, XS_file=XS_file, masks_file=masks_file)
        xs_df = get_cs_df_from_gdf(self.reach)
        self.update_cs_gdf(xs_df, XS_file)
        return xs_df

    def PointToMesh(self,
                    num_vertices: int|None = None,
                    space_min: float|None = None,
                    space_max: float|None = None, 
                    with_original: bool = False, 
                    output_point_file: str|None = None,
                    output_line_file: str|None = None,  
                    bins: int|None = None,
                    
                    positions: list[float]|None = None, 
                    bw: list[float]|float|None = None, 
                    dist_threshold: float = 2.0, 
                    width_threshold: float =0.0, 
                    extend_to_boundary: bool = True,

                    output_dir: str|None = None, 
                    XS_file: str|None = None, 
                    masks_file: str|None = None, 
                    extend_boundary: bool = False, 
                    buffer: float = 0, 
                    input_DEM_file: str|None = None, 
                    output_DEM_file: str|None = None

                    ) -> pd.DataFrame:

        """This function is used to generate mesh files (points or line strings) from a observed points file.

        Parameters
        ----------
        num_vertices : int, optional
            The number of vertices in each cross-section, by default None
        space_min : float, optional
            The min spacing between cross-sections, by default None
        space_max : float, optional
            The max spacing between cross-sections, by default None
        with_original : bool, optional
            To extract values from observed points after interpolating cross-sections or not, by default False
        output_point_file : str, optional
            The path to save the mesh files as line features. by default None
        output_line_file : str, optional
            The path to save the mesh files as linestring features. by default None
        bins : int, optional
            The number of bins for the signal function to find out the locations of cross-sections, by default None
        extend_to_boundary : bool, optional
            To extend the interpolation range to the boundary or not, by default True
        buffer : int, optional
            The size of buffer to extract the values from the original observed points, by default 0
        positions : list[float], optional
            Given the positions of cross-sections by S-coordinates, by default None
        bw : list[float]/float, optional
            The width of buffers that could be set as one value (float) or as different values (list), by default None
        dist_threshold : float, optional
            The maximum distance between two cross-sections.
            If two cross-sections are larger than this threshold, a cross-section will be inserted to the middel, by default 2.0
        width_threshold : float, optional
            The ratio of the cross-section length and channel width should larger than this threshold, or it will be removed, by default 0.0
        output_dir : str, optional
            The path of directory to store the output files, by default None
        XS_file : str, optional
            The file name or path of the cross-section file.
            If output_dir is specified, it will be saved in output_dir, by default None
        masks_file : str, optional
            The file name or path of the cross-section file.
            If output_dir is specified, it will be saved in output_dir, by default None
        input_DEM_file : str, optional
            The file path of the DEM file for integrating the mesh with DEM, by default None
        output_DEM_file : str, optional
            The file path of the raster file that is integration of mesh and DEM, by default None

        Returns
        -------
        DataFrame
            The DataFrame stores the x, y, S, N, and z coordinates of vertices in each cross-section.
            The XS_id in the table indicates which cross-section each point belongs to.
        """
  
        output_point_file = file_mangement.get_output_path(output_dir, output_point_file)
        output_line_file = file_mangement.get_output_path(output_dir, output_line_file)

        xs_df = self.PointToXS(bins, num_vertices, positions, bw, dist_threshold=dist_threshold, width_threshold=width_threshold, extend_to_boundary=extend_to_boundary, output_dir=output_dir, XS_file=XS_file, masks_file=masks_file)
        allvertices = curvilinear_interp.generate_mesh(reach=self.reach,
                                                       num_vertices=num_vertices,
                                                       xs_df=xs_df, space_min=space_min,
                                                       space_max=space_max,
                                                       with_original=with_original,
                                                       output_point_file=output_point_file,
                                                       output_line_file=output_line_file,
                                                       extend_boundary=extend_boundary,
                                                       buffer=buffer,
                                                       input_DEM_file=input_DEM_file,
                                                       output_DEM_file=output_DEM_file)
        return allvertices

    def __interpXS(self, xs_up, xs_dw, pos, i, num_vertices):
        width_up = xs_up['N'].max() - xs_up['N'].min()
        # get the normal vector
        xs_temp = np.array(self.reach.coord_converter.sn2xy_coord(pos, 1, 0))
        xs_mid_p = np.array(self.reach.coord_converter.sn2xy_coord(pos, 0, 0))
        N_vec = xs_temp - xs_mid_p
        # extend line along with the normal direction
        xs_bd_p = [xs_mid_p + N_vec * w for w in (-width_up, width_up)]
        xs_line = LineString(xs_bd_p)
        xs_line = self.reach.bd_geom.intersection(xs_line)
        # get N coordinates of line clipped by the boundary
        xs_bd_N = [np.linalg.norm(xs_mid_p - p) for p in xs_line.coords]
        xs_bd_N[0] *= -1

        S_up = xs_up['S'].unique()[0]
        S_dw = xs_dw['S'].unique()[0]
        prop = (S_dw - pos) / (S_dw - S_up)

        xs_up.reset_index(drop=True, inplace=True)
        xs_dw.reset_index(drop=True, inplace=True)
        xs_md = pd.DataFrame({'XS_id':i, 'S': pos, 'N':np.linspace(*xs_bd_N, num_vertices)})
        xs_md['z_up'] = xs_up['z']
        xs_md['z_dw'] = xs_dw['z']
        xs_md['z'] = xs_md['z_up'] * prop + xs_md['z_dw'] * (1 - prop)
        xs_md.drop(columns=['z_up', 'z_dw'], inplace=True)
        xs_md['z'] = xs_md['z'].ffill()
        xs_md['z'] = xs_md['z'].bfill()
        return xs_md

    def MeshToXS(self,
                 positions: list[float]|None = None,
                 space: float|None = None,
                 num_vertices: int|None = None,
                 interp: str|None = 'linear',
                 output_dir: str|None = None,
                 output_file: str|None = None
                 ):
        """This function is used to interpolate cross-sections at given locations from a given mesh.

        Parameters
        ----------
        positions : list[float], required
            The S coordinates of cross-section locations, by default None
        space : float, optional
            A spacing to create cross-sections uniformly, by default None
        num_vertices : int, optional
            The number of vertices in each cross-section, by default None
        interp : str, optional
            The method to do interpolation based on the mesh, by default 'linear'
        output_dir : str, optional
            The path of directory to store the output files, by default None
        output_file : str, optional
            The name the output file (path: output_dir/output_file), by default None

        Returns
        -------
        DataFrame
            The DataFrame stores the x, y, S, N, and z coordinates of vertices in each cross-section.
            The XS_id in the table indicates which cross-section each point belongs to.
        """
        
        print('Interpolating cross-sections from mesh data...')
        output_dir = output_dir if output_dir else self.reach.output_dir
        file_mangement.check_dir_exists(output_dir)
        output_file = file_mangement.get_output_path(output_dir, output_file)

        ## preprocessing input data
        ## input data may be mesh points or lines
        xs_df = self.reach.ms_gdf
        geom_type = xs_df.geometry[0].geom_type
        
        def getPositions(xs_df, positions = None, space = None):
            '''
            generating cross-sections at
                the positions given by user or,
                the positions by a even space from start to end
            '''
            S_ls = xs_df['S'].unique()
            if positions is None:
                if space:
                    min_S, max_S = min(S_ls), max(S_ls)
                    positions = np.arange(min_S, max_S, space)
                else:
                    raise ValueError("Please give positions or a space for getting cross-sections!")    
            return positions

        if (geom_type == 'Point'):  
            ## set the number of vertices: 1. from user 2. from original mesh point
            xs_id_ls = xs_df['XS_id'].unique()
            num_vertices = num_vertices if num_vertices else xs_df[xs_df['XS_id'] == xs_id_ls[0]]['N'].count()
            xs_df = self.reach.convert_sn_coord_for_layer(xs_df)
            g = xs_df.groupby(by='XS_id')
            new_S = g['S'].mean()
            xs_df.drop(columns='S', inplace=True)
            xs_df = xs_df.join(new_S, on='XS_id', how='inner')
            positions = getPositions(xs_df, positions, space)
            S_ls = xs_df['S'].unique()
    
            dfs = []
            for i, pos in enumerate(positions):
                if pos in S_ls:
                    xs_md = xs_df[xs_df['S'] == pos]
                else:
                    idx = (np.abs(S_ls - pos)).argmin()
                    idx = max(0, idx-1)
                    S_up, S_dw = S_ls[idx:idx+2]
                    xs_up = xs_df[xs_df['S'] == S_up]
                    xs_dw = xs_df[xs_df['S'] == S_dw]
                    xs_md = self.__interpXS(xs_up, xs_dw, pos, i, num_vertices)
                dfs.append(xs_md)

        elif (geom_type == 'LineString'):
            ## filter out the lines which are parallel to the centerline
            ## the cross-sections should be perpendicular to the centerline
            xs_df = xs_df[xs_df.geometry.intersection(self.reach.cl_geom).geom_type == 'Point']
            ints_bd_l = xs_df.geometry.intersection(self.reach.bd_geom).length
            max_l = ints_bd_l.max()
            median_l = ints_bd_l.median()
            cl_l = self.reach.cl_geom.length
            while( abs(max_l - cl_l) < abs(max_l - median_l)):
                xs_df = xs_df[ints_bd_l != max_l]
                ints_bd_l = xs_df.geometry.intersection(self.reach.bd_geom).length
                max_l = ints_bd_l.max()
                median_l = ints_bd_l.median()
                
            ## get corresponding S coord to each cross-section
            _ = get_cs_df_from_gdf(self.reach, xs_df) # S coord will be added to the input df by this function 
            positions = getPositions(xs_df, positions, space)      
            S_min = xs_df['S'].unique().min()
            S_max = xs_df['S'].unique().max()

            dfs = []
            for i, pos in enumerate(positions):
                if ((pos < S_min) or (pos > S_max)):
                    continue
                p = self.reach.coord_converter.sn2xy_point(s=pos, n=0, z=0)
                xs_df['dist'] = xs_df.geometry.distance(p)
                xs_df.sort_values(by='dist', inplace=True)
                xs_df.reset_index(inplace=True)
                xs_geom = xs_df.loc[:1, ['index', 'geometry']]
                xs_up_geom = xs_geom[xs_geom['index'] == xs_geom['index'].min()]['geometry'].values[0]
                xs_dw_geom = xs_geom[xs_geom['index'] == xs_geom['index'].max()]['geometry'].values[0]
                xs_df.drop(columns='index', inplace=True)
                ## use getXSDataFrame??
                def getDFfromLineGeom(line):
                    '''
                        get XYZ coords and SNZ coords from LineString geometry
                        to do interpolation
                    '''
                    xyz_array = np.array(line.coords)
                    sn_array = preprocess.convert_xy2sn(xyz_array[:, :2], self.cl_xyz, self.bw_org)
                    coords_array = np.concatenate([xyz_array, sn_array], axis=1)
                    columns = {'x', 'y', 'z', 'S', 'N'}
                    temp_df = pd.DataFrame(coords_array, columns=columns)
                    return temp_df 
                
                xs_up = getDFfromLineGeom(xs_up_geom)
                xs_dw = getDFfromLineGeom(xs_dw_geom)
                num_vertices = num_vertices if num_vertices else xs_up['N'].count()
                xs_md = self.__interpXS(xs_up, xs_dw, pos, i, num_vertices)
                dfs.append(xs_md)
                
        ## combining all cross-sections together
        xs_df = pd.concat(dfs, ignore_index=True)
        if self.reach.elev_base :
            xs_df['z'] = self.reach.elev_base - xs_df['z']
        if self.reach.shift:
            self.loelev_baset_elev = xs_df['z'].min()
            xs_df['z'] = xs_df['z'] - self.loelev_baset_elev
        getXYbyRow = lambda row: self.reach.coord_converter.sn2xy_point(*row)
        xs_df['geometry'] = xs_df[['S', 'N', 'z']].apply(getXYbyRow, axis=1)
        xs_df = gpd.GeoDataFrame(xs_df, geometry='geometry', crs=self.reach.proj_crs)
        file_mangement.export_geofile(xs_df, output_file)
        return xs_df
