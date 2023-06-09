"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import scipy.sparse as sp
import numpy as np
from collections import Counter

from rph_kmeans.point_reducer_base import RPPointReducerBase
from rph_kmeans._point_reducer_cy import get_ary_labels
from rph_kmeans._point_reducer_cy import update_densex_and_weight, update_sparsex_and_weight, update_labels
from rph_kmeans._point_reducer_cy import densex_radius_bkt_improve, sparsex_radius_bkt_improve
from rph_kmeans.utils import cal_dist2_ary_sparse, cal_dist2_ary_dense


class RPPointReducerCy(RPPointReducerBase):
	def __init__(self, w=None, max_point=2000, proj_num=5, max_iter=1000, sample_dist_num=1000,
			bkt_improve=None, radius_divide=None, bkt_size_keepr=1.0, center_dist_keepr=1.0, verbose=1):
		super(RPPointReducerCy, self).__init__(w, max_point, proj_num, max_iter, sample_dist_num,
			bkt_improve, radius_divide, bkt_size_keepr, center_dist_keepr, verbose=verbose)


	def fit_transform(self, X):
		"""
		Args:
			X (numpy.ndarray or scipy.sparse.csr_matrix): (sample_num, feature_size)
		Returns:
			np.ndarray: (reduced_point_num, feature_size); reduced points
			np.ndarray: (reduced_point_num,); weight, number of samples belonging to each reduced point
			np.ndarray: (sample_num,); labels, indicating which reduced point each sample belongs to
			int: number of iteration
		"""
		X = self.check_input_X(X)
		self.check_max_point(self.max_point, X)
		self.cal_dist2_ary_func = cal_dist2_ary_sparse if self.sparse_x else cal_dist2_ary_dense

		sample_num, feature_size = X.shape
		self.get_w(X)
		iter_num = 0
		group_num = sample_num

		labels = np.arange(0, sample_num, dtype=np.uint32)
		group_weight = np.ones((sample_num,), dtype=X.dtype)
		reduced_X = X

		if self.verbose > 0:
			print('Iteration begin: X.shape = {}, max_point = {}, w = {}, proj_num = {}'.format(X.shape, self.max_point, self.w, self.proj_num))
		while iter_num < self.max_iter:
			if group_num <= self.max_point:
				if self.verbose > 0:
					print('Reduced point number {} <= max_point {}. Iteration stop.'.format(group_num, self.max_point))
				break

			proj_vecs, b = self.gen_proj(feature_size, self.w)
			pj_mat = self.random_projection(reduced_X, proj_vecs, b)  # (group_num, projection_num)
			bkt_ary, bkt_num = get_ary_labels(pj_mat)

			if self.bkt_improve == 'radius':
				bkt_num = self.radius_bkt_improve(reduced_X, bkt_ary, self.radius2, bkt_num)
			elif self.bkt_improve == 'min_bkt_size':
				if self.bkt_size_keepr < 1.0:
					bkt_num = self.min_bkt_size_bkt_improve(bkt_ary)
			elif self.bkt_improve == 'min_center_dist':
				if self.center_dist_keepr < 1.0:
					bkt_num = self.min_center_dist_bkt_improve(bkt_ary, reduced_X)
			else:
				if self.bkt_improve is not None:
					raise RuntimeError("Parameter bkt_improve must be one of {None, 'radius', 'min_bkt_size', 'min_center_dist'}")

			update_labels(labels, bkt_ary)
			reduced_X, group_weight = self.update_x_and_weight(reduced_X, group_weight, bkt_ary, bkt_num)

			group_num = bkt_num
			iter_num += 1

			if iter_num == 1:
				self.cal_dist2_ary_func = cal_dist2_ary_dense

			if self.verbose > 1:
				group_bkt_count, orphan_point_count = self.count_group_orphan(bkt_ary)
				print('Iter {}: Reduced point number = {}; Group bucket number={}; Orphan point number={}'.format(
					iter_num, group_num, group_bkt_count, orphan_point_count))

		if self.verbose > 0:
			print('Total iteration = {}; Number of reduced points = {}'.format(iter_num, group_num))
		return reduced_X, group_weight.flatten(), labels, iter_num


	def radius_bkt_improve(self, reduced_X, bkt_ary, R2, bkt_num):
		if sp.issparse(reduced_X):
			data, indices, indptr = reduced_X.data, reduced_X.indices, reduced_X.indptr
			return sparsex_radius_bkt_improve(data, indices, indptr, reduced_X.shape[0], reduced_X.shape[1], bkt_ary, R2, bkt_num)
		else:
			return densex_radius_bkt_improve(reduced_X, bkt_ary, R2, bkt_num)


	def bkt_ary_buckets_keepr_wrapper(self, bkt_ary, handle_func, *args, **kwargs):
		buckets = self.bkt_ary_to_buckets(bkt_ary)
		group_buckets, orphan_points = self.split_group_orphan(buckets)

		sorted_idx, keep_num = handle_func(group_buckets, *args, **kwargs)

		for i, idx in enumerate(sorted_idx[:keep_num]):
			bkt_ary[group_buckets[idx]] = i
		orphan_points.extend([p_idx for bkt_idx in sorted_idx[keep_num:] for p_idx in group_buckets[bkt_idx]])
		bkt_num = len(orphan_points) + keep_num
		bkt_ary[orphan_points] = np.arange(keep_num, bkt_num, dtype=np.uint32)

		return bkt_num


	def min_bkt_size_bkt_improve(self, bkt_ary):
		# TODO: write with c++
		return self.bkt_ary_buckets_keepr_wrapper(bkt_ary, self.min_bkt_size_bkt_improve_)


	def min_bkt_size_bkt_improve_(self, group_buckets):
		sorted_idx = np.argsort([len(bkt) for bkt in group_buckets])
		keep_num = int(len(group_buckets) * self.bkt_size_keepr)
		return sorted_idx, keep_num


	def min_center_dist_bkt_improve(self, bkt_ary, X):
		# TODO: write with c++
		return self.bkt_ary_buckets_keepr_wrapper(bkt_ary, self.min_center_dist_bkt_improve_, X)


	def min_center_dist_bkt_improve_(self, group_buckets, X):
		keep_num = int(len(group_buckets) * self.center_dist_keepr)
		dist_ary = []
		for bkt in group_buckets:
			center = np.mean(X[bkt], axis=0)
			if sp.issparse(X):
				center = sp.csr_matrix(center)
			dist_ary.append(np.median(self.cal_dist2_ary_func(X[bkt], center)))  # (bucket_size,)
		sorted_idx = np.argsort(dist_ary)
		return sorted_idx, keep_num


	def update_x_and_weight(self, reduced_X, weight, bkt_ary, bkt_num):
		if sp.issparse(reduced_X):
			return update_sparsex_and_weight(reduced_X, weight, bkt_ary, bkt_num)
		else:
			return update_densex_and_weight(reduced_X, weight, bkt_ary, bkt_num)


	def bkt_ary_to_buckets(self, bkt_ary):
		d = {}
		for i, bkt_id in enumerate(bkt_ary):
			d.setdefault(bkt_id, []).append(i)
		return list(d.values())


	def count_group_orphan(self, bkt_ary):
		counter = Counter(bkt_ary)
		orphan_count = 0
		for bkt_id, member_count in counter.items():
			if member_count == 1:
				orphan_count += 1
		return len(counter) - orphan_count, orphan_count


