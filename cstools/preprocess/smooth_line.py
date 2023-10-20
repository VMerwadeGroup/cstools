from __future__ import annotations

import numpy as np
from scipy import interpolate
from shapely.geometry import LineString

def smooth_line_geom(line_geom: LineString, spacing: float|None = 10) -> LineString:
    num = int(line_geom.length / spacing) + 1
    xy = np.array(line_geom.coords)
    xy = xy[:, :2]
    tck, u = interpolate.splprep(xy.T, u=None, k=3, s=0) 
    u_new = np.linspace(u.min(), u.max(), num)
    new_x, new_y = interpolate.splev(u_new, tck)
    return LineString(zip(new_x, new_y))

    
