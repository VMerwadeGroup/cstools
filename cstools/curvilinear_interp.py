from __future__ import annotations

from .reach import *
from .file_mangement.export import export_geofile
from shapely import speedups
from shapely.ops import nearest_points

speedups.enable()

def _line_spanning_bd(p1: Point, p2: Point, reach: RiverReach) -> LineString:
    """Find a line including p1 and p2 spans the boundary

    Parameters
    ----------
    p1 : Point
        A point close to the left bank.
    p2 : Point
        A point close to the right bank.
    reach : RiverReach
        RiverReach object provides the centerline geometry, boundary geometry, and the approximated channel width.

    Returns
    -------
    LineString
        A Line spans the boundary.
    """
    p_on_cl = nearest_points(p1, reach.cl_geom)[1]
    p_on_cl_a = np.array(p_on_cl.xy)
    p1_a = np.array(p1.xy)
    p2_a = np.array(p2.xy)
    
    ## p1 and p2 is used to provide the unit norm vector
    v = (p2_a - p1_a)
    v /= np.linalg.norm(v, ord=2)

    ## make rays start from the point on the centerline and point to left/right normal direction
    ## the length of the rays are set as 5 times reach.approx_width
    left_ext = LineString([p_on_cl_a, p_on_cl_a + 5*reach.approx_width*v])
    right_ext= LineString([p_on_cl_a, p_on_cl_a - 5*reach.approx_width*v])

    ## the end points are nearests points from the intersections of the boundary to the point on the centerline.
    left_inters = left_ext.intersection(reach.bd_geom.boundary)
    left_end = nearest_points(left_inters, p_on_cl)[0]

    right_inters = right_ext.intersection(reach.bd_geom.boundary)
    right_end = nearest_points(right_inters, p_on_cl)[0]
    return LineString([left_end, right_end])        
    

def process_csdf_for_mesh(reach: RiverReach,
                          xs_df: pd.DataFrame,
                          num_vertices: float,
                          extend_boundary: bool = False,
                          buffer: float = 0,
                          input_DEM_file: str|None = None,
                          output_DEM_file: str|None = None,
                          output_dir: str = './temp'
                          ) -> pd.DataFrame:

    ## divide cross-sections in to the same proportions
    indice = xs_df['XS_id'].unique()
    group = xs_df.groupby('XS_id')
    starts = group['N'].min()
    ends = group['N'].max()
    locations = group['S'].mean()

    ### TODO: make the cross-section can access values from DEM raster

    # if (extend_boundary):
    #     buffer = self.approx_width * 0.25 if (buffer == 0) else buffer
    #     if (input_DEM_file):
    #         DEM_raster = rioxarray.open_rasterio(input_DEM_file)
    #     else:
    #         DEM_raster = self.getDEMRaster(output_dir, output_DEM_file)
    #     self.bd_geom = self.bd_geom.buffer(buffer)
    #     self.bd_gdf.loc[0, 'geometry'] = self.bd_geom

    print('Generating points on given corss-sections...')
    vertices_dfs =[]
    for i, s, start_N, end_N in zip(indice, locations, starts, ends):
        filt = xs_df[xs_df['XS_id'] == i]
        filt = filt.sort_values(by=['N'])
        filt.reset_index(drop=True, inplace=True)

        if reach.bd_geom:
            # find the N of intersections between cross-sections and the boundary
            p1 = Point(filt.loc[0, ['x','y']].values)
            p2 = Point(filt.loc[len(filt)-1, ['x','y']].values)
            sp_line = _line_spanning_bd(p1, p2, reach)
            bd_SN = [reach.coord_converter.xy2sn_coord(*c) for c in sp_line.coords]
            bd_N = sorted([SN[1] for SN in bd_SN])
        else:
            bd_N = sorted([start_N, end_N])

        if (bd_N[0]*bd_N[1] > 0):
            print(i, bd_N, filt)
    
        vertices_N = np.linspace(*bd_N, num_vertices)
        f = interpolate.interp1d(filt['N'], filt['z'], kind='linear', bounds_error=False, fill_value=np.nan)
        vertices_z = f(vertices_N)
        dic = {'XS_id': i, 'S':s, 'N':vertices_N, 'z':vertices_z, 'vertex_id':range(num_vertices)}
        vertices_df = pd.DataFrame(dic)
        
        ### TODO: make the cross-section can access values from DEM raster
        # getXY = lambda row: self.coord_converter.sn2xy_coord(*row)

        # if (extend_boundary):
        #     xyz = vertices_df[['S', 'N']].apply(getXY, axis=1)
        #     x, y, _ = zip(*xyz)
        #     vertices_df['x'] = x
        #     vertices_df['y'] = y

        #     points_df = self.getRasterValue(points=filt, DEM_raster=DEM_raster)
        #     if self.shift:
        #         points_df['raster_value'] = points_df['raster_value']/self.vertical_conversion_factor - self.lowest_elev
        #     vertices_df['z_DEM'] = points_df['raster_value']
        #     if (np.isnan(vertices_df.loc[0, 'z'])):
        #         vertices_df.loc[0, 'z'] = vertices_df.loc[0, 'z_DEM']
        #     if (np.isnan(vertices_df.loc[len(vertices_df)-1, 'z'])):
        #         vertices_df.loc[len(vertices_df)-1, 'z'] = vertices_df.loc[len(vertices_df)-1, 'z_DEM']
        #     vertices_df.drop(columns=['x', 'y', 'z_DEM'], inplace=True)

        vertices_df['z'] = vertices_df['z'].interpolate(method='linear')
        vertices_df['z'] = vertices_df['z'].fillna(method='bfill')
        vertices_df['z'] = vertices_df['z'].fillna(method='ffill')
        vertices_dfs.append(vertices_df)

    return pd.concat(vertices_dfs)

def fix_mesh_intersection(reach: RiverReach, mesh_points: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    g = mesh_points.groupby(['XS_id'])
    mesh_S = pd.DataFrame(g['geometry'].apply(lambda x: LineString(x.tolist())))
    mesh_S['S'] = g['S'].median()
    mesh_S.sort_values(by=['S'], inplace=True)
    xs_df = gpd.GeoDataFrame(mesh_S, crs=reach.proj_crs)
    xs_df = xs_df.reset_index()
    xs_df['FID'] = xs_df.index

    num_inter = []
    inter_idx = []
    for i in range(len(xs_df)):
        xs = xs_df.iloc[i]['geometry']
        inter = xs_df[xs_df['FID'] != i]['geometry'].intersects(xs)
        num_inter.append(sum(inter))
        inter[i] = False
        inter_idx.append(xs_df[inter].index)
        
    xs_df['num_inter'] = num_inter
    xs_df['inter_idx'] = inter_idx
    xs_inter = xs_df[xs_df['num_inter'] != 0]

    while (len(xs_inter) > 0):
        row = xs_inter.iloc[0]
        up_idx = row['FID']
        visited = []
        inter_ls = list(row['inter_idx'])
        while (len(inter_ls) > 0):
            inter_idx = inter_ls.pop()
            if (inter_idx not in visited):
                visited.append(inter_idx)
                next_xs = xs_inter.loc[inter_idx]
                inter_ls += list(next_xs['inter_idx'])
        xs_inter = xs_inter.drop(visited)
        dw_idx = max(visited)

        up_idx = up_idx - 1 if up_idx > 0 else up_idx
        dw_idx = dw_idx + 1 if dw_idx < xs_df.index[-1] else xs_df.index[-1]
        xs_up = xs_df.iloc[up_idx]
        xs_dw = xs_df.iloc[dw_idx]
        bd = reach.bd_geom.boundary
        dists = [bd.project(Point(xs.geometry.coords[idx])) for xs in [xs_up, xs_dw] for idx in [0, -1]]
        
        z_up = np.array([p[2] for p in xs_up.geometry.coords])
        z_dw = np.array([p[2] for p in xs_dw.geometry.coords])
        
        num_xs = dw_idx - up_idx + 1
        num_gaps = num_xs - 1
        num_vertices = len(z_up)

        left = np.linspace(dists[0], dists[2], num_xs)[1:-1]
        right = np.linspace(dists[1], dists[3], num_xs)[1:-1]
        # xs_df.iloc[up_idx:dw_idx+1].plot()

        for i, (l, r) in enumerate(zip(left, right)):
            idx = up_idx + i + 1
            s = (xs_up['S'] * (num_gaps-(i+1)) + xs_dw['S'] * (i+1)) / num_gaps

            line = LineString([bd.interpolate(l), bd.interpolate(r)])
            positions = np.linspace(0, line.length, num_vertices)
            points = [line.interpolate(pos) for pos in positions]
            z = tuple((z_up * (num_gaps-(i+1)) + z_dw * (i+1)) / num_gaps)
            xyz = [(*p.coords[0], z) for p,z  in zip(points, z)]
            x, y ,z = zip(*xyz)
            new_line = LineString(xyz)
            xs_df.loc[idx, 'S'] = s
            xs_df.loc[idx, 'geometry'] = new_line

            ## modify the mesh points
            snz = [reach.coord_converter.xy2sn_coord(xi, yi) for xi, yi in zip(x,y)]
            s, n, _ = zip(*snz)
            xs_id = xs_df.loc[idx, 'XS_id']
            filter_idx = mesh_points[mesh_points['XS_id'] == xs_id].index

            dic = {
                'XS_id': xs_id,
                'geometry': [Point(p) for p in xyz],
                'x': x,
                'y': y,
                'z': z,
                'S': s,
                'N': n
            }

            df = pd.DataFrame(dic, index=filter_idx)
            mesh_points.update(df, join='left', overwrite=True)
            
        # xs_df.iloc[up_idx:dw_idx+1].plot()
    mesh_S = xs_df[['XS_id', 'geometry']]
    return mesh_S, mesh_points

def generate_mesh(reach: RiverReach,
                 xs_df: pd.DataFrame,
                 num_vertices: float|None = None,
                 space_min: float|None = None,
                 space_max: float|None = None,
                 equal_space: bool = False,
                 keep_last: bool = False,
                 with_original: bool = False,
                 output_point_file: str|None = None,
                 output_line_file: str|None = None,
                 extend_boundary: bool = False,
                 buffer: float = 0,
                 input_DEM_file: str|None = None,
                 output_DEM_file: str|None = None
                 ) -> gpd.GeoDataFrame:
    '''
    generating mesh from cross-sections by interpolating through SN coordinate
    Input:
        xs_df: dataframe of cross-sections
        num_vertices: number of vertices for each cross-section
        space_min: min spacing betwenn cross-sections
        space_max: max spacing betwenn cross-sections
        equal_space: spacing between cross-sections should be equal or not
        keep_last: if the last interpolated cross-section is less than space_min, keep it or not
        with_original: output data should include input cross-sections or only interpolated cross-se ctions
        output_file: provide path to output shapefile of interpolated cross-sections
    Output:
        A dataframe of interpolated cross-sections
        A shapefile of interpolated cross-sections
    '''
    
    ## check the values of interpolation spaces
    if ((space_min is None) and (space_max is None)):
        space_min = reach.interp_xs_space
        space_max = reach.interp_xs_space
    elif ((space_min is None) and (space_max is not None)):
        space_min = space_max
    elif ((space_min is not None) and (space_max is None)):
        space_max = space_min
    
    if (num_vertices is None):
        if (extend_boundary):
            if (buffer == 0):
                num_vertices = int(reach.vertices_in_xs * 1.5)
            else:
                num_vertices = int((reach.approx_width + buffer * 2) / reach.bw_org)
        else:
            num_vertices = reach.vertices_in_xs

    dfInverseSN = lambda x: reach.coord_converter.sn2xy_point(*x)
    
    proportion = np.linspace(0, 1, num_vertices)
    group = xs_df.groupby('XS_id')
    locations = group['S'].mean()
    s_diff = locations.diff()

    given_vertices = process_csdf_for_mesh(reach=reach, xs_df=xs_df, num_vertices=num_vertices,
                                              extend_boundary=extend_boundary, buffer=buffer,
                                              input_DEM_file=input_DEM_file, output_DEM_file=output_DEM_file)
    nums = np.ceil(s_diff/space_max)-1 # num + 1 = # of segments, num = # of lines

    interp_dfs = []
    print('Interpolating cross-sections...')
    for idx in locations.index[1:]:
    # for idx in tqdm(locations.index[1:], desc='Interpolating cross-sections...'):
        up_s, down_s = locations[idx-1], locations[idx]
        num = int(nums[idx])

        # detemine the S coordinate of vertices
        if equal_space:
            num += keep_last - 1
            if ((up_s+space_max) - (down_s-space_max)) < space_min:
                new_s = np.linspace(up_s+space_max, down_s-space_max, num)
            else:
                new_s = np.linspace(up_s+space_max, down_s-space_max, num)
        elif (down_s - (up_s + num * space_max)) >= space_min:
            new_s = np.linspace(up_s+space_max, up_s+space_max*num, num)
        else:
            num += keep_last - 1
            new_s = np.linspace(up_s+space_max, up_s+space_max*num, max(num,0))

        for i, s in enumerate(new_s):
            upstream = given_vertices[given_vertices['XS_id'] == idx-1]
            downstream = given_vertices[given_vertices['XS_id'] == idx]
            vertices = upstream.join(downstream, on='vertex_id', how='right',\
                lsuffix='_up', rsuffix='_down')
                
            ## z: linear interpolation between upstream and downstream cross-sections
            dist_up = np.abs(s - vertices['S_up'].unique()[0])
            dist_down = np.abs(s - vertices['S_down'].unique()[0])
            weight_up = dist_down / (dist_up + dist_down)
            weight_down = dist_up / (dist_up + dist_down)
            vertices['z'] = vertices['z_up'] * weight_up + vertices['z_down'] * weight_down

            ## s: given
            vertices['S'] = s

            ## n: spanning line
            p1 = reach.coord_converter.sn2xy_point(s=s, n=-1)
            p2 = reach.coord_converter.sn2xy_point(s=s, n=1)

            spanning_line = _line_spanning_bd(p1, p2, reach)
            start_p, end_p = spanning_line.coords
            start_SN = reach.coord_converter.xy2sn_coord(*start_p)
            end_SN = reach.coord_converter.xy2sn_coord(*end_p)

            bd_n = sorted([start_SN[1], end_SN[1]])
            n_arr = np.interp(proportion, [0, 1], bd_n)
            vertices['N'] = n_arr

            ## setup a dataframe
            interp_df = vertices.loc[:, ['S', 'N', 'z', 'vertex_id']]
            interp_df['geometry'] = interp_df[['S', 'N', 'z']].apply(dfInverseSN, axis=1)
            interp_df['XS_id'] = f'{idx-1}_{i}'

            ## get values from the original observed points
            if with_original:
                ob_gdf = reach.ob_gdf
                ob_gdf.sort_values(by=['S', 'N'], inplace=True)
                ob_gdf.reset_index(drop=True, inplace=True)

                mask = spanning_line.buffer(reach.bw_org)
                selected = ob_gdf[ob_gdf.geometry.within(mask)]
                selected.dropna(inplace=True)

                def comInfoInOb(row, selected):
                    if (not selected.empty):
                        selected['distance'] = selected.geometry.distance(row.geometry)
                        selected.sort_values('distance', inplace=True, ascending=True)
                        if len(selected) >= 12:
                            selected = selected.iloc[:12]
                        if selected['distance'].min() < 1e-8:
                            z = selected[selected['distance'] == selected['distance'].min()]['z']
                        else:
                            selected['weight'] = 1 / selected['distance'] ** 2
                            z = (selected['z'] * selected['weight']).sum() / selected['weight'].sum()
                    else:
                        z = row['z']
                    return z
                interp_df['z'] = interp_df.apply(comInfoInOb, axis=1, args=(selected,))
            interp_dfs.append(interp_df)
    

    if len(interp_dfs) > 0:
        interp_vertices = pd.concat(interp_dfs)
        interp_vertices.sort_values(by=['S', 'N'], inplace=True)
        interp_vertices.reset_index(drop=True, inplace=True)
        given_vertices['geometry'] = given_vertices[['S', 'N', 'z']].apply(dfInverseSN, axis=1)
        allvertices = pd.concat([given_vertices, interp_vertices])
    else:
        given_vertices['geometry'] = given_vertices[['S', 'N', 'z']].apply(dfInverseSN, axis=1)
        allvertices = given_vertices
    allvertices['XS_id'] = allvertices['XS_id'].astype(str)
    allvertices = allvertices.sort_values(by=['XS_id', 'vertex_id'])

    ## mesh points
    mesh_points = gpd.GeoDataFrame(allvertices, geometry='geometry')                
    mesh_points.crs = reach.proj_crs
    mesh_points.sort_values(by=['S', 'N'], inplace=True)
    mesh_points.reset_index(drop=True, inplace=True)
    mesh_line_S, mesh_points = fix_mesh_intersection(reach, mesh_points)

    ## mesh lines
    # mesh_S = mesh_points.groupby(['XS_id'])['geometry'].apply(lambda x: LineString(x.tolist()))
    # mesh_line_S = gpd.GeoDataFrame(mesh_S, crs=self.proj_crs) 
    mesh_N = mesh_points.groupby(['vertex_id'])['geometry'].apply(lambda x: LineString(x.tolist()))
    mesh_line_N = gpd.GeoDataFrame(mesh_N, crs=reach.proj_crs) 

    mesh_lines = pd.concat([mesh_line_S, mesh_line_N])
    mesh_lines.reset_index(drop=True, inplace=True)

    if output_point_file:
        export_geofile(mesh_points, output_point_file)
    if output_line_file:
        export_geofile(mesh_lines, output_line_file)
        
    return mesh_points
