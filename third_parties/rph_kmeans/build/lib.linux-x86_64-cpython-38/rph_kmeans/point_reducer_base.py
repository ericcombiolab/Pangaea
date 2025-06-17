"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import warnings
import numpy as np
import scipy.sparse as sp
from sklearn.metrics.pairwise import paired_distances
from sklearn.utils import check_array

from rph_kmeans.utils import check_return


class RPPointReducerBase(object):
	def __init__(self, w=None, max_point=2000, proj_num=5, max_iter=1000, sample_dist_num=1000,
			bkt_improve=None, radius_divide=None, bkt_size_keepr=1.0, center_dist_keepr=1.0, verbose=False):
		self.w = w
		self.max_point = max_point
		self.proj_num = proj_num
		self.max_iter = max_iter
		self.sample_dist_num = int(sample_dist_num)
		self.bkt_improve = bkt_improve
		self.radius2 = radius_divide**2 if radius_divide is not None else None
		self.bkt_size_keepr=bkt_size_keepr
		self.center_dist_keepr = center_dist_keepr
		self.verbose=verbose
		self.sparse_x = None


	def fit_transform(self, X):
		raise NotImplementedError


	@check_return('sample_dist')
	def get_sample_dist(self, X):
		return paired_distances(
			X[np.random.choice(X.shape[0], self.sample_dist_num)],
			X[np.random.choice(X.shape[0], self.sample_dist_num)],
			'euclidean'
		)


	@check_return('w')
	def get_w(self, X):
		sample_dist = self.get_sample_dist(X)
		w = np.median(sample_dist) * 0.5
		print(f'Note: w of RPH is automatically set to {w}')
		return w


	def gen_proj(self, dim, w):
		b = np.random.uniform(0, 1, size=(self.proj_num,))
		proj_vecs = np.random.normal(0, 1.0/w, size=(dim, self.proj_num))  # (feature_size, projection_num)
		return proj_vecs, b


	def random_projection(self, features, proj_vecs, b):
		return (features.dot(proj_vecs) + b).astype(np.int32) # (sample_num, projection_num)


	def check_input_X(self, X):
		X = check_array(X, accept_sparse="csr", order='C', dtype=[np.float64, np.float32])
		self.sparse_x = sp.issparse(X)
		return X


	def check_max_point(self, max_point, X):
		if max_point >= X.shape[0]:
			warnings.warn("max_point is larger than sample number of input X. The process of point reduction won't run")


	def split_group_orphan(self, buckets):
		"""
		Returns:
			list: group_buckets; [[point_idx, ...], ...]
			list: orphan_points; [point_idx, ...]
		"""
		group_buckets, orphan_points = [], []
		for bkt in buckets:
			if len(bkt) == 1:
				orphan_points.append(bkt[0])
			else:
				group_buckets.append(bkt)
		return group_buckets, orphan_points

