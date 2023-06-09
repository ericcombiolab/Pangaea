"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import scipy.sparse as sp
import numpy as np

from rph_kmeans.point_reducer_base import RPPointReducerBase
from rph_kmeans.utils import cal_dist2_ary_sparse, cal_dist2_ary_dense
from rph_kmeans.utils import cal_weighted_ave_sparse_vec, cal_weighted_ave_dense_vec

class RPPointReducerPy(RPPointReducerBase):
	def __init__(self, w=None, max_point=2000, proj_num=5, max_iter=1000, sample_dist_num=1000,
			bkt_improve=None, radius_divide=None, bkt_size_keepr=1.0, center_dist_keepr=1.0, verbose=1):
		super(RPPointReducerPy, self).__init__(w, max_point, proj_num, max_iter, sample_dist_num,
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
		self.cal_weighted_ave_vec_func = cal_weighted_ave_sparse_vec if self.sparse_x else cal_weighted_ave_dense_vec
		self.reduced_op_after_func = lambda X: X.A1 if self.sparse_x else lambda X: X

		sample_num, feature_size = X.shape
		self.get_w(X)
		iter_num = 0
		group_num = sample_num

		labels = np.arange(0, sample_num, dtype=np.uint32)
		group_weight = np.ones((sample_num, 1), dtype=X.dtype)
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
			buckets = self.pjmat_to_buckets(pj_mat)  # [[point_idx, ...], ...];
			group_buckets, orphan_points = self.split_group_orphan(buckets)

			if self.bkt_improve == 'radius':
				assert self.radius2 is not None
				group_buckets, orphan_points = self.radius_bkt_improve(group_buckets, orphan_points, reduced_X)
			elif self.bkt_improve == 'min_bkt_size':
				if self.bkt_size_keepr < 1.0:
					group_buckets, orphan_points = self.min_bkt_size_bkt_improve(group_buckets, orphan_points)
			elif self.bkt_improve == 'min_center_dist':
				if self.center_dist_keepr < 1.0:
					group_buckets, orphan_points = self.min_center_dist_bkt_improve(group_buckets, orphan_points, reduced_X)
			elif self.bkt_improve is None:
				pass
			else:
				raise RuntimeError(
					"Parameter bkt_improve = {}. It must be one of "
					"[None, 'radius', 'min_bkt_size', 'min_center_dist']".format(self.bkt_improve))

			group_num = len(group_buckets) + len(orphan_points)
			labels = self.update_labels(labels, reduced_X, group_buckets, orphan_points, group_num)
			reduced_X, group_weight = self.update_x_and_weight(reduced_X, group_weight, group_buckets, orphan_points, group_num)

			iter_num += 1

			if iter_num == 1:
				self.cal_dist2_ary_func = cal_dist2_ary_dense
				self.cal_weighted_ave_vec_func = cal_weighted_ave_dense_vec

			if self.verbose > 1:
				print('Iter {}: Reduced point number = {}; Group bucket number={}; Orphan point number={}'.format(
					iter_num, group_num, len(group_buckets), len(orphan_points)))
		if self.verbose > 0:
			print('Total iteration = {}; Number of reduced points = {}'.format(iter_num, group_num))
		return reduced_X, group_weight.flatten(), labels, iter_num


	def pjmat_to_buckets(self, pj_mat):
		"""
		Args:
			pj_mat (np.ndarray): (sample_num, projection_num)
		Returns:
			list: [bucket1, bucket2, ...], bucket=[point_idx, ...]
		"""
		d = {}
		for i in range(pj_mat.shape[0]):
			d.setdefault(pj_mat[i].tobytes(), []).append(i)
		return d.values()


	def radius_bkt_improve(self, group_buckets, orphan_points, X):
		new_group_buckets = []
		for bkt in group_buckets:
			_, groups = self.radius_single_bkt_improve(X[bkt], self.radius2, self.cal_dist2_ary_func)
			for g in groups:
				if len(g) == 1:
					orphan_points.append(bkt[g[0]])
				else:
					new_group_buckets.append([bkt[r] for r in g])
		return new_group_buckets, orphan_points


	def radius_single_bkt_improve(self, X, R2, dist2_ary_func):
		centers = [0]
		groups = [[0]]
		for i in range(1, X.shape[0]):
			dist2_ary = dist2_ary_func(X[centers], X[i])
			c = np.argmin(dist2_ary)
			if dist2_ary[c] < R2:
				groups[c].append(i)
			else:
				centers.append(i); groups.append([i])
		return centers, groups


	def min_bkt_size_bkt_improve(self, group_buckets, orphan_points):
		sorted_idx = np.argsort([len(bkt) for bkt in group_buckets])
		keep_num = int(len(group_buckets) * self.bkt_size_keepr)
		new_group_buckets = [group_buckets[idx] for idx in sorted_idx[:keep_num]]
		orphan_points.extend([p_idx for bkt_idx in sorted_idx[keep_num:] for p_idx in group_buckets[bkt_idx]])
		return new_group_buckets, orphan_points


	def min_center_dist_bkt_improve(self, group_buckets, orphan_points, X):
		keep_num = int(len(group_buckets) * self.center_dist_keepr)
		dist_ary = []
		for bkt in group_buckets:
			center = X[bkt].mean(axis=0)
			if sp.issparse(X):
				center = sp.csr_matrix(center)
			dist_ary.append(np.median(self.cal_dist2_ary_func(X[bkt], center)))  # (bucket_size,)
		sorted_idx = np.argsort(dist_ary)
		new_group_buckets = [group_buckets[idx] for idx in sorted_idx[:keep_num]]
		orphan_points.extend([p_idx for bkt_idx in sorted_idx[keep_num:] for p_idx in group_buckets[bkt_idx]])
		return new_group_buckets, orphan_points


	def update_labels(self, labels, reduced_X, group_buckets, orphan_points, new_group_num):
		labels_bridge = np.zeros(reduced_X.shape[0], dtype=np.uint32)
		for i, bkt in enumerate(group_buckets):
			labels_bridge[bkt] = i
		labels_bridge[orphan_points] = np.arange(len(group_buckets), new_group_num)
		labels = labels_bridge[labels]
		return labels


	def update_x_and_weight(self, reduced_X, group_weight, group_buckets, orphan_points, new_group_num):
		shrink_X_new, group_weight_new = [], np.zeros((new_group_num, 1), dtype=reduced_X.dtype)
		for i, bucket_ranks in enumerate(group_buckets):
			weight = group_weight[bucket_ranks]
			group_weight_new[i] = np.sum(weight)
			weight /= group_weight_new[i]
			ave_vec = self.cal_weighted_ave_vec_func(reduced_X[bucket_ranks], weight)
			shrink_X_new.append(ave_vec)
		shrink_X_new.append(reduced_X[orphan_points].toarray() if sp.issparse(reduced_X) else reduced_X[orphan_points])
		group_weight_new[len(group_buckets): new_group_num] = group_weight[orphan_points]
		reduced_X, group_weight = np.vstack(shrink_X_new), group_weight_new
		return reduced_X, group_weight



