from __future__ import annotations

from .reach import *
from . import preprocess
from . import file_mangement

def get_cs_df_from_gdf(reach: RiverReach, cs_file: str|None = None) -> pd.DataFrame:
    '''
    get cooridnates and attributes from layer
    make data as pandas dataframe
    '''
    if isinstance(cs_file, gpd.GeoDataFrame):
        gdf = cs_file
    elif type(cs_file) == str:
        gdf = file_mangement.init_cross_sections(cs_file)
        reach.cs_gdf = gdf
    elif cs_file is not None:
        raise TypeError('''cs_file should be a path (str) or a GeoDataFrame''') 
    elif reach.cs_gdf is not None:
        gdf = reach.cs_gdf
    else:
        raise ValueError('''Please run getXSbyHist() to get XS3D file, 
        or input a GeoDataFrame or a vector shapefile layer''') 

    gdf['geometry'] = gdf['geometry'].intersection(reach.bd_geom)
    gdf.reset_index(drop=True, inplace=True)
    gdf['intersection'] = gdf['geometry'].intersection(reach.cl_geom)
    gdf = gdf[~gdf['intersection'].is_empty]
    ### Note: it could be multiple intersections instead of just one
    
    SNZ = gdf['intersection'].apply(lambda ints: reach.coord_converter.xy2sn_coord(*(ints.coords[0])))
    S, _, _ = zip(*SNZ)
    gdf['S'] = S
    gdf.reset_index(drop=True, inplace=True)
    dfs = []

    ## explode each line into points
    for i, row in gdf.iterrows():
        points = row.geometry.coords
        SN = [reach.coord_converter.xy2sn_coord(*p) for p in points]
        (_, n, z) = zip(*SN)
        x, y, _ = zip(*points)
        dic = {'XS_id':i, 'x': x, 'y':y, 'S':row['S'], 'N': n, 'z':z}
        df = pd.DataFrame(dic)
        dfs.append(df)
    
    xs_df = pd.concat(dfs, ignore_index=True)
    xs_df.sort_values(by=['S', 'N'])
    if reach.ob_gdf is None:
        if reach.elev_base:
            xs_df['z'] = reach.elev_base - xs_df['z']
        if reach.shift:
            reach.lowest_elev = xs_df['z'].min()
            xs_df['z'] = xs_df['z'] - reach.lowest_elev
    
    xs_df['z'] = xs_df['z'].apply(lambda z: np.nan if (np.isnan(z) or z < -1e5) else z)
    if xs_df['z'].isnull().sum() > 0:
        group = xs_df.groupby('XS_id')
        xs_df['z'] = group['z'].apply(lambda g: g.interpolate(method='linear')).reset_index()['z']
        xs_df['z'] = group['z'].ffill().reset_index()['z']
        xs_df['z'] = group['z'].bfill().reset_index()['z']
    return xs_df

def project_linear_survey(row):
    cs = row['geometry']
    clip_ob = row['clip_obs']

    # projecting observation points to the segments of a cross-section
    N = clip_ob.apply(cs.project)
    points = N.apply(cs.interpolate)
    xy = points.apply(lambda row: row.coords[0])
    z = clip_ob.apply(lambda row: row.z)
    x, y = zip(*xy.values)
    df = pd.DataFrame(zip(x, y, z, N), columns=['x', 'y', 'z', 'N'])
    df.sort_values(by='N', inplace=True)
    line = LineString(df[['x', 'y', 'z']].to_numpy())
    row['Points'] = len(clip_ob)
    row['geometry'] = line
    return row

def project_zigzag_survey(row, num):
    # get the start points and end points in cross-section shapefile
    cs = row['geometry']
    clip_ob = row['clip_obs']
    length = cs.length
    vertices_N = np.linspace(0, length, num) 
    bin_cover = np.concatenate(([0], np.convolve(vertices_N, [0.5, 0.5], 'valid'), [length]))

    # projecting observation points to the segments of a cross-section
    df = pd.DataFrame()
    df['N'] = clip_ob.apply(cs.project)
    df['distance'] = clip_ob.apply(cs.distance)
    df['z'] = clip_ob.apply(lambda row: row.z)
                
    start_point, end_point = cs.coords[0], cs.coords[-1]
    seg = LineString([start_point, end_point])

    new_points = []
    cnt_z_nan = 0
    for j in range(len(bin_cover)-1):
        filt = (df['N'] > bin_cover[j]) & (df['N'] <= bin_cover[j+1]) 
        selected = df[filt]
        selected.dropna(inplace=True)
        selected.reset_index(drop=True, inplace=True)
        if not selected.empty:
            selected.sort_values('distance', inplace=True, ascending=True)
            if len(selected) >= 12:
                selected = selected.iloc[:12]
            if selected['distance'].min() < 1e-8:
                z = selected[selected['distance'] == selected['distance'].min()]['z']
            else:
                selected['weight'] = 1 / selected['distance'] ** 2
                z = (selected['z'] * selected['weight']).sum() / selected['weight'].sum()
        else:
            z = -1e39
            cnt_z_nan += 1
        new_p = seg.interpolate(vertices_N[j])
        new_points.append((*new_p.coords[0], z))

    ## make sure there is at least one value in a cross-section line
    if cnt_z_nan != num:
        line = LineString(new_points)
        row['geometry'] = line
    else:
        row['geometry'] = np.nan
    return row

def project_points_cs(reach: RiverReach,
                      bws: list[float],
                      num: float|None = None
                      ) -> None:
    '''
    project observations to cross-section seperately

    given layers of observations and cross-sections
    return cross-sections with elevations in a layer
    '''
    print('Projecting observation points to cross-sections lines...')

    # create point objects that are projection on the staright cross-section
    
    cs_gdf = reach.cs_gdf.to_crs(reach.proj_crs)
    ob_gdf = reach.ob_gdf.to_crs(reach.proj_crs)
    
    masks = cs_gdf.buffer(np.array(list(bws)))
    cs_gdf['clip_obs'] = [ob_gdf[ob_gdf.geometry.within(mask)]['geometry'] for mask in masks.geometry]

    if reach.survey_type == 'linear':
        cs_gdf = cs_gdf.apply(project_linear_survey, axis=1)
    elif reach.survey_type == 'zigzag':
        num = num if num else reach.vertices_in_xs
        cs_gdf = cs_gdf.apply(project_zigzag_survey, args=(num,), axis=1)
    cs_gdf = cs_gdf[~cs_gdf['geometry'].isnull()]
    cs_gdf = cs_gdf.drop(columns='clip_obs')

    reach.cs_gdf= cs_gdf
    reach.masks_gdf = masks

def check_cs_distance(reach: RiverReach,
                      positions: np.array,
                      widths: np.array,
                      dist_threshold: float,
                      ) -> tuple[np.array, np.array]:
    '''
    A private function used for finding positions by histogram.
    This function is going to make sure the distance between 2 CS is not too large.
    If the distance is larger than the threshold, then add a new CS in the middle.

    postions, widths: from the __findPositions function based on histogram

    threshold: dist_threshold (defalut value = 2.0) * channel width
    new width: mean of the original widths
    '''

    area = reach.bd_geom.area
    cl_geom_in_bd = reach.bd_geom.intersection(reach.cl_geom)
    length = cl_geom_in_bd.length
    threshold = dist_threshold * area / length
    avg_width = np.mean(widths)

    positions = np.array(positions)
    pos_diff = np.abs(np.diff(positions))
    idx_over = np.argwhere(pos_diff > threshold)
    idx_over = np.reshape(idx_over, (-1,))
    while (len(idx_over) > 0):
        add_pos = [np.mean(positions[idx:idx+2]) for idx in idx_over]
        positions = np.insert(positions, idx_over+1, add_pos)
        widths = np.insert(widths, idx_over+1, avg_width)

        pos_diff = np.abs(np.diff(positions))
        idx_over = np.argwhere(pos_diff > threshold)
        idx_over = np.reshape(idx_over, (-1,))
    return positions, widths

def find_positions_by_hist(reach: RiverReach,
                           series,
                           bins,
                           dist_threshold,
                           extend_to_boundary = True
                           ) -> dict[float, float]:
    '''
    find out peaks of observations histogram
    use np.histogram to generate histogram arrays
    use scipy.signal to find out peaks in the arrays
    return the locations and widths of peaks

    Input: pandas series of observations (S,N), bins given by user
    Output: locations of peaks, widths of peaks 
    '''
    print('Determine positions of cross-sections by the histogram of observation points...')
    division = np.linspace(series.min(), series.max(), bins)
    counts, division = np.histogram(series, bins=division)
    
    ## scipy.signal.find_peaks will not recognize a peak if it is located at the end 
    counts_bf = np.concatenate(([counts[1]], counts , [counts[-2]]))
    criteria = np.mean(counts_bf)/np.sqrt(bins)
    peaks, prop = signal.find_peaks(counts_bf, prominence=criteria)
    
    ## peak_widths gives the widths of indice to counts, 
    widths, heights, left, right = signal.peak_widths(counts_bf, peaks, rel_height=0.60)
    ## transfer indice widths into division widths
    widths *= division[1] - division[0]

    ## index should -1 due to adding boundary,
    ## but +0.5 due to indice in counts is 0.5 large than in division
    left -= 0.5
    right -= 0.5
    
    ## interp_to_S = lambda x: np.interp(x, range(len(division)), division)
    left_bd = np.interp(left, range(len(division)), division)
    right_bd = np.interp(right, range(len(division)), division)

    ## positions = np.convolve(division, [0.5, 0.5], mode='valid')
    centroids = [series[(series <= r) & (series >= l)].mean() for l, r in zip(left_bd, right_bd)]
    peaks -= 1 # move -1 bakc for the original series
    widths = np.where(widths == 0, np.mean(widths), widths)

    ## the cross-sectional survey (line does not need to check the distance between cross-sections)
    if (reach.survey_type == 'zigzag'):
        if (extend_to_boundary):
            centroids = [series.min()] + centroids + [series.max()]
            avg_width = np.mean(widths)
            widths = np.insert([avg_width, avg_width], [1], widths)
        if (dist_threshold):
            centroids, widths = check_cs_distance(reach, centroids, widths, dist_threshold)

    return dict(zip(centroids, widths))

def locate_cs_by_hist(reach = RiverReach,
                      bins = None,
                      points_in_cs = None,
                      positions = None,
                      bw = None, 
                      extend_to_boundary = True, 
                      dist_threshold=2., 
                      width_threshold=0., 
                      XS_file = None, 
                      masks_file = None):
    '''
    Determine the better locations of cross-sections and projected observations on to them
    Input:
        bins: the number of bins that is used in searching positions of cross-sections 
        points_in_cs: the number of points that is included in each cross-section line (only for zigzag survey)
        bw: buffer widths for each mask that is used to clip observations for projecting
    Output: cross-sections(LineStringZ) and masks(Polygon) vector layer
    '''
    ## check the values of parameters
    if bins is None:
        bins = round(np.sqrt(len(reach.ob_gdf['S'])))

    if points_in_cs is None:
        points_in_cs = reach.vertices_in_xs

    if positions is None:
        cs_positions_dict = find_positions_by_hist(reach, reach.ob_gdf['S'], bins, dist_threshold, extend_to_boundary)
        if len(cs_positions_dict) < 4:
            space = preprocess.get_avg_channel_width(reach.bd_geom, reach.cl_geom)
            S_max, S_min = reach.ob_gdf['S'].max(), reach.ob_gdf['S'].min()
            length = S_max - S_min
            if (length < space * 10):
                space = length / 10
            centroids = np.arange(S_min + 0.02 * length, S_max - 0.02 * length, space)
            bw = space / 4
        else:
            centroids = cs_positions_dict.keys()

        if bw:
            widths = bw if hasattr(bw, '__iter__') else [bw] * len(centroids)
            cs_positions_dict = dict(zip(centroids, widths))
    else:
        widths = bw if hasattr(bw,  '__iter__') else [bw] * len(positions)
        cs_positions_dict = dict(zip(positions, widths))
        
    reach.bw_median = np.median(list(cs_positions_dict.values()))
    ## create a cross-section layer
    print('Creating cross-section lines...')

    cross_sections = []
    bws = []
    for (p, b) in cs_positions_dict.items():
        filt = (reach.ob_gdf['S'] <= p + b) & (reach.ob_gdf['S'] >= p - b)
        selected = reach.ob_gdf[filt]
        n_max, n_min = selected['N'].max(), selected['N'].min()
        
        ## check the cross-section exists
        if (np.isnan(n_max) or np.isnan(n_min)):
            continue

        ## check the cross-section is long enough
        n_large = np.max(np.abs([n_max, n_min]))
        cs = LineString([reach.coord_converter.sn2xy_point(s=p, n=n) for n in [n_min, n_max]])
        extended_cs = LineString([reach.coord_converter.sn2xy_point(s=p, n=n*10) for n in [-n_large, n_large]])
        intersect_cs = extended_cs.intersection(reach.bd_geom)
        ch_width = intersect_cs.length
        cs_width = cs.length

        if (cs_width > width_threshold * ch_width):
            cross_sections.append(cs) 
            bws.append(b)

    reach.cs_gdf = gpd.GeoDataFrame(geometry=cross_sections, crs=reach.proj_crs)
    project_points_cs(reach, bws, points_in_cs)
    reach.cs_gdf.geometry = reach.cs_gdf.geometry.set_crs(reach.proj_crs)
    
    # XS_path = os.path.join(self.output_dir, XS_file)
    if XS_file:
        file_mangement.export_geofile(reach.cs_gdf, XS_file)

    if masks_file:
        file_mangement.export_geofile(reach.masks_gdf, masks_file)