"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import warnings
import numpy as np
from copy import deepcopy
from sklearn.cluster import KMeans

class RPHKMeans(object):
	def __init__(self, n_clusters=8, n_init=1, point_reducer_version='cy',
			w=None, max_point=2000, proj_num=5, max_iter=1000, sample_dist_num=1000,
			bkt_improve=None, radius_divide=None, bkt_size_keepr=1.0, center_dist_keepr=1.0,
			reduced_kmeans_kwargs=None, final_kmeans_kwargs=None, verbose=1):
		"""
		Args:
			n_clusters (int): default: 8
				The number of clusters to form as well as the number of centroids to generate.
			n_init (int): default: 1
				Number of time the rp k-means algorithm will be run. The final results will be
				the best output of n_init consecutive runs in terms of inertia.
			point_reducer_version (str): {'cy', 'py'}, default: 'cy'
				Version of point reducer module. 'cy' is for cython version and 'py' is for python version.
				If cython version is failed to import, python version will be used instead.
			w (float): default: None
				Width of bucket for random projection process. If set to None, w will be the half of
				the median of distances of paired samples from input X.
			max_point (int): default: 2000
				Algorithm will stop when reduced point number is less than max_point. If max_point is larger than
				or equal to sample number of X, the process of point reduction won't run.
			proj_num (int): default: 5
				Number of vector for random projection process.
			max_iter (int): default: 1000
				Maximum number of iterations of the algorithm to run.
			sample_dist_num (int): default: 1000
				Number of paired samples chosen to decide default w. It will be ignored when w is set by user.
			bkt_improve (str or None): {None, 'radius', 'min_bkt_size', 'min_center_dist'}
				Methods of improving bucket quality.
			radius_divide (float): default: None
				Radius for 'radius' bucket-improving method. If set to None, no bucket-improving process
				will be ran even though bkt_improve is set to 'radius'. See paper or code for details.
			bkt_size_keepr (float): default: 0.8
				Keep ratio for 'min_bkt_size' bucket-improving method. If set to 1.0, no bucket-improving
				process will be ran even though bkt_improve is set to 'min_bkt_size'. See paper or code for details.
			center_dist_keepr (float): default: 0.8
				Keep ratio for 'min_center_dist' bucket-improving method. If set to 1.0, no bucket-improving
				process will be ran even though bkt_improve is set to 'min_center_dist'. See paper or code for details.
			reduced_kmeans_kwargs (dict): default: None
				kwargs of kmeans to find centers of reduced point. If set to None, default kwargs will be used.
			final_kmeans_kwargs (dict): default: None
				kwargs of kmeans after center initialization. If set to None, default kwargs will be used.
			verbose (int): {0, 1, 2}, default: 1
				Controls the verbosity. Print nothing when set to 0 and print most details when set to 2.

		Attributes:
			cluster_centers_ (np.ndarray): (n_clusters, n_features)
				Coordinates of cluster centers.
			labels_ (np.ndarray): (n_samples)
				Labels of each point
			inertia_ (float):
				Sum of squared distances of samples to their closest cluster center.
			n_iter_ (int):
				Number of kmeans's iterations run after center initialization.

			reduced_X_ (np.ndarray): (n_reduced_point, n_features)
				Reduced points.
			reduced_X_weight_ (np.ndarray): (n_reduced_point,)
				Weight of reduced points. reduced_X_weight_[i] represents the number of original points merged into
				reduced_X_[i].
			rp_labels_ (np.ndarray): (n_samples,)
				Labels of each original point indicating which reduced point it belongs to.
				Specifically， X[i] belongs to reduced_X_[rp_labels_[i]].
			init_centers_ (np.ndarray): (n_clusters, n_features)
				Initial centers.
			rp_iter_ (int):
				Number of iteration of point reduction process.
		"""
		if point_reducer_version == 'cy':
			try:
				from rph_kmeans.point_reducer_cy import RPPointReducerCy
				RPPointReducer = RPPointReducerCy
			except:
				warnings.warn('The cython version of rph-kmeans is not installed properly. Use python version instead.')
				from rph_kmeans.point_reducer_py import RPPointReducerPy
				RPPointReducer = RPPointReducerPy
		else:
			from rph_kmeans.point_reducer_py import RPPointReducerPy
			RPPointReducer = RPPointReducerPy
		self.point_reducer = RPPointReducer(w, max_point, proj_num, max_iter, sample_dist_num,
			bkt_improve, radius_divide, bkt_size_keepr, center_dist_keepr, verbose)
		self.n_clusters = n_clusters
		self.n_init = n_init
		self.reduced_kmeans_kwargs = deepcopy(reduced_kmeans_kwargs) if reduced_kmeans_kwargs is not None else {}
		if 'n_clusters' not in self.reduced_kmeans_kwargs:
			self.reduced_kmeans_kwargs['n_clusters'] = n_clusters
		self.final_kmeans_kwargs = deepcopy(final_kmeans_kwargs) if final_kmeans_kwargs is not None else {}
		if 'n_clusters' not in self.final_kmeans_kwargs:
			self.final_kmeans_kwargs['n_clusters'] = n_clusters
		self.final_kmeans_kwargs['n_init'] = 1

		self.reduced_X_ = None
		self.reduced_X_weight_ = None
		self.rp_labels_ = None
		self.rp_iter_ = None
		self.init_centers_ = None
		self.rp_inertia_ = None

		self.final_kmeans_clt = None
		self.inertia_ = None
		self.cluster_centers_ = None
		self.labels_ = None
		self.n_iter_ = None


	def init_centers(self, X):
		reduced_X, group_weight, labels, rp_iter = self.point_reducer.fit_transform(X)
		assert len(reduced_X) == len(group_weight)
		if len(reduced_X) < self.n_clusters:
			raise RuntimeError('Number of reduced points is too small, please try smaller w or larger proj_num')
		reduced_clt = KMeans(**self.reduced_kmeans_kwargs)
		while True:
			try:
				y_pred = reduced_clt.fit_predict(reduced_X, sample_weight=group_weight)
			except IndexError:
				print('Index error raised, re-run kmeans for skeleton to initialize centers. This may due to bugs in sklearn.')
			else:
				break
		return reduced_clt.cluster_centers_, reduced_X, group_weight, labels, rp_iter, reduced_clt.inertia_, y_pred


	def fit_predict(self, X):
		"""
		Args:
			X (numpy.ndarray or scipy.sparse.csr_matrix): (n_samples, n_features)
		Returns:
			np.ndarray: (n_samples,)
		"""
		self.fit(X)
		return self.labels_


	def fit(self, X):
		"""
		Args:
			X (numpy.ndarray or scipy.sparse.csr_matrix): (n_samples, n_features)
			Training instances to cluster. It must be noted that the data will
			be converted to C ordering, which will cause a memory copy
			if the given data is not C-contiguous.
		"""
		self.inertia_ = np.inf
		for i in range(self.n_init):
			init_centers_, reduced_X_, reduced_X_weight_, rp_labels_, rp_iter_, rp_inertia_, rp_y_pred_ = self.init_centers(X)
			clt = KMeans(init=init_centers_, **self.final_kmeans_kwargs)
			clt.fit(X)
			if clt.inertia_ < self.inertia_:
				self.inertia_ = clt.inertia_
				self.final_kmeans_clt = clt
				self.init_centers_, self.reduced_X_, self.reduced_X_weight_, self.rp_labels_, self.rp_iter_, self.rp_inertia_, self.rp_y_pred_ = \
					init_centers_, reduced_X_, reduced_X_weight_, rp_labels_, rp_iter_, rp_inertia_, rp_y_pred_
		self.cluster_centers_ = self.final_kmeans_clt.cluster_centers_
		self.labels_, self.n_iter_ = self.final_kmeans_clt.labels_, self.final_kmeans_clt.n_iter_


	def predict(self, X):
		"""
		Args:
			X (numpy.ndarray or scipy.sparse.csr_matrix): (n_samples, n_features)
		Returns:
			np.ndarray: (n_samples,)
		"""
		X = check_array(X, accept_sparse="csr", order='C', dtype=[np.float64, np.float32])
		return self.final_kmeans_clt.predict(X)
