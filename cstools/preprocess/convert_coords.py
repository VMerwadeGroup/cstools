from __future__ import annotations

import numpy as np
from shapely.geometry import Point, LineString
import geopandas as gpd
# import jax.numpy as jnp
# from jax import jit

def get_xy_array(ob_gdf: gpd.GeoSeries) -> np.array:
    getXYZ = lambda geom: geom.coords[0]
    XYZ = ob_gdf.apply(getXYZ)
    x, y, _ = zip(*XYZ)
    return np.array([x, y]).T

def get_xyz_array(ob_gdf: gpd.GeoSeries) -> np.array:
    getXYZ = lambda geom: geom.coords[0]
    XYZ = ob_gdf.apply(getXYZ)
    x, y, z = zip(*XYZ)
    return np.array([x, y, z]).T

def get_xy_array_from_line(l_geom: LineString) -> np.array:
	l_array = np.array(l_geom.coords)
	l_array = l_array[:, :2] # remove "z" coordinate, if it has
	return l_array

def get_xyz_array_from_line(l_geom: LineString) -> np.array:
	l_array = np.array(l_geom.coords)
	return l_array

def convert_xy2sn(
		xy: np.array, 
		cl_array: np.array, 
		bw_org: float
	) -> np.array:
	
	num_cl_p = cl_array.shape[0]
	xy_stack = np.repeat(xy[:, np.newaxis, :], num_cl_p-1, axis=1)

	cl_vec = np.diff(cl_array, axis=0)
	cl_len = np.linalg.norm(cl_vec, ord=2, axis=1)
	avg_seg_len = np.mean(cl_len)

	point_start_vec = (xy_stack - cl_array[:num_cl_p-1])
	point_start_dist = np.linalg.norm(point_start_vec, ord=2, axis=2)
	area = np.cross(point_start_vec, cl_vec)

	h = area / cl_len
	h_abs = np.abs(h)

	t_stack = np.sum(np.multiply(point_start_vec, cl_vec), axis=2) / np.square(cl_len)
	filt = np.logical_and(np.logical_and(t_stack<=1.00, t_stack>=-0.00), h_abs < (avg_seg_len + bw_org))
	h_abs_filt = filt * h_abs
	row_sum = np.sum(h_abs_filt, axis=1)
	acute_row = np.where(~(row_sum == 0))
	obtuse_row = np.where((row_sum == 0))

	obtuse_dist = point_start_dist[obtuse_row]
	p_idx = np.nanargmin(obtuse_dist, axis=1)
	cl_cum_len = np.append(np.array([0]), np.cumsum(cl_len))
	s_obtuse = cl_cum_len[p_idx]
	sign = np.sign(area[obtuse_row][range(p_idx.shape[0]), p_idx])
	n_obtuse = sign * np.min(obtuse_dist, axis=1)

	h_abs_acute = h_abs_filt[acute_row]
	t_stack_acute = t_stack[acute_row]
	h_abs_acute[np.isclose(h_abs_acute,0)] = np.nan
	idx_acute = np.nanargmin(h_abs_acute, axis=1)
	t_pick = t_stack_acute[np.arange(idx_acute.shape[0]),idx_acute]

	seg_len = cl_len[idx_acute]
	s_start = cl_cum_len[idx_acute]
	s_acute = s_start + (seg_len * t_pick)
	n_acute = h[acute_row][np.arange(idx_acute.shape[0]),idx_acute]

	sn = np.zeros_like(xy) 
	sn[acute_row,0] = s_acute
	sn[acute_row,1] = n_acute
	sn[obtuse_row,0] = s_obtuse
	sn[obtuse_row,1] = n_obtuse
	
	return sn

def convert_sn2xy(
		sn: np.array, 
		cl_array: np.array, 
		bw_org: float
	) -> np.array:

	num_cl_p = cl_array.shape[0]
	cl_vec = np.diff(cl_array, axis=0)
	cl_len = np.linalg.norm(cl_vec, ord=2, axis=1)
	cl_len_cum = np.cumsum(cl_len)
	
	

# ---
# JAX implementations

# @jit
# def get_t_jax(point_start_vec, cl_vec):
#     cl_len = jnp.linalg.norm(cl_vec, ord=2, axis=1)
#     t_stack = jnp.sum(jnp.multiply(point_start_vec, cl_vec), axis=2) / jnp.square(cl_len)
#     return t_stack

# def convert_xy2sn_jax(
# 		xy: np.array, 
# 		cl_array: np.array, 
# 		bw_org: float
# 	) -> np.array:
	
# 	num_cl_p = jnp.size(cl_array,0)
# 	xy_stack = jnp.repeat(xy[:, jnp.newaxis, :], num_cl_p-1, axis=1)

# 	cl_vec = jnp.diff(cl_array, axis=0)
# 	cl_len = jnp.linalg.norm(cl_vec, ord=2, axis=1)
# 	avg_seg_len = jnp.mean(cl_len)

# 	point_start_vec = (xy_stack - cl_array[:num_cl_p-1])
# 	point_start_dist = jnp.linalg.norm(point_start_vec, ord=2, axis=2)
# 	area = jnp.cross(point_start_vec, cl_vec)

# 	h = area / cl_len
# 	h_abs = jnp.abs(h)
		
# 	t_stack = get_t_jax(point_start_vec, cl_vec)
# 	filt = jnp.logical_and(jnp.logical_and(t_stack<=1.00, t_stack>=-0.00), h_abs < (avg_seg_len + bw_org))
# 	h_abs_filt = filt * h_abs
# 	row_sum = jnp.sum(h_abs_filt, axis=1)
# 	acute_row = jnp.nonzero(row_sum != 0)
# 	obtuse_row = jnp.nonzero(row_sum == 0)

# 	obtuse_dist = point_start_dist[obtuse_row]
# 	idx_obtuse = jnp.nanargmin(obtuse_dist, axis=1).astype(int)
# 	cl_cum_len = jnp.append(jnp.array([0]), jnp.cumsum(cl_len))
# 	s_obtuse = cl_cum_len[idx_obtuse]
# 	sign = jnp.sign(area[obtuse_row][jnp.arange(idx_obtuse.shape[0]), idx_obtuse])
# 	n_obtuse = sign * jnp.min(obtuse_dist, axis=1)

# 	h_abs_acute = h_abs_filt[acute_row]
# 	t_stack_acute = t_stack[acute_row]
# 	h_abs_acute = jnp.where(jnp.isclose(h_abs_acute,0), jnp.nan, h_abs_acute)
# 	idx_acute = jnp.nanargmin(h_abs_acute, axis=1).astype(int)
# 	t_pick = t_stack_acute[jnp.arange(idx_acute.shape[0]),idx_acute]

# 	seg_len = cl_len[idx_acute]
# 	s_start = cl_cum_len[idx_acute]
# 	s_acute = s_start + (seg_len * t_pick)
# 	n_acute = h[acute_row][jnp.arange(idx_acute.shape[0]),idx_acute]

# 	sn = np.zeros_like(xy) 
# 	sn[acute_row,0] = s_acute
# 	sn[acute_row,1] = n_acute
# 	sn[obtuse_row,0] = s_obtuse
# 	sn[obtuse_row,1] = n_obtuse

# 	return sn

# ---
# geom function for iteration

class CoordConverter(object):
	def __init__(self, cl_geom: LineString):
		self.cl_geom = cl_geom
		self.line_verts = [Point(p) for p in cl_geom.coords]
		self.line_dist = np.array([cl_geom.project(p) for p in self.line_verts])
	
	def xy2sn_coord(self,
					x: float|None = None, 
					y: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		"""Covert (x, y, z) to (s, n, z) with a centerline
		The point will be used is a Point object is given.
		Otherwise, (x, y, z) should be specified.

		Parameters
		----------
		cl_geom : LineString
			The geometry of the centerline (flowline).
		x : float | None, optional
			x coordinate of the point, by default None
		y : float | None, optional
			y coordinate of the point, by default None
		z : float | None, optional
			z coordinate of the point, by default None
		point : Point | None, optional
			Point object with (x, y, z), by default None

		Returns
		-------
		tuple[float, float, float|None]
			3D: (s, n, z) or 2D: (s, n, None)
		"""

		if (point == None):
			point = Point(x, y, z) if z else Point(x, y)
		s = self.cl_geom.project(point)
		n = self.cl_geom.distance(point)
		z = point.z if point.has_z else None

		if (s > self.cl_geom.length):
			intp_prev, intp_later = self.line_verts[-2:]
		elif (s <= 0):
			intp_prev, intp_later = self.line_verts[0:2]
		else:
			idx_ls = self.line_dist[(self.line_dist - s) < 0]
			start_idx = (s - idx_ls).argmin()
			intp_prev, intp_later = self.line_verts[start_idx:start_idx+2]

		# # Ring: previous point -> interpolated point -> point
		# lr = LinearRing([intp_prev, intp_later, point])       
		
		# # S-N: cw = right N+, ccw = left N-
		# if lr.is_ccw:
		#     n *= -1

		flow_vec = np.array([intp_later.x - intp_prev.x, intp_later.y - intp_prev.y])
		point_vec = np.array([point.x - intp_prev.x, point.y - intp_prev.y])
		n *= np.sign(np.cross(point_vec, flow_vec)) # sign of z value of cross-product

		return s, n, z

	def xy2sn_point(self,
					x: float|None = None, 
					y: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		s, n, z = self.xy2sn_coord(x, y, z, point)
		point = Point(s, n) if z is None else Point(s, n, z)
		return point


	def sn2xy_coord(self,
					s: float|None = None, 
					n: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		"""_summary_

		Parameters
		----------
		cl_geom : LineString
			The geometry of the centerline (flowline).
		s : float | None, optional
			s coordinate of the point, by default None
		n : float | None, optional
			n coordinate of the point, by default None
		z : float | None, optional
			z coordinate of the point, by default None
		point : Point | None, optional
			Point object with (s, n, z), by default None

		Returns
		-------
		tuple[float, float, float|None]
			3D: (x, y, z) or 2D: (x, y, None)
		"""
		if point:
			s, n, z = point.coords[0]
		intp = self.cl_geom.interpolate(s)
		if (s > self.cl_geom.length):
			intp_prev = self.line_verts[-1]
		elif (s < 0):
			intp_prev = self.line_verts[0]
		elif (s == 0):
			intp_prev, intp = self.line_verts[0], self.line_verts[1]
		else:
			idx_ls = self.line_dist[(self.line_dist - s) < 0]
			start_idx = (s - idx_ls).argmin()
			intp_prev = self.line_verts[start_idx]
		dist = intp.distance(intp_prev)
		# get unit normal vector of the segment
		vec_N = np.array([intp.y - intp_prev.y, intp_prev.x - intp.x]) / dist
		x, y = intp.coords[0] + vec_N * n
		return x, y, z

	def sn2xy_point(self,
					s: float|None = None, 
					n: float|None = None,
					z: float|None = None,
					point: Point|None = None
					) -> tuple[float, float, float|None]:
		x, y, z = self.sn2xy_coord(s, n, z, point)
		point = Point(x, y) if z is None else Point(x, y, z)
		return point