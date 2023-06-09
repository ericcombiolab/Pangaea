"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import scipy
from kneed import KneeLocator
from multiprocessing import Pool, cpu_count
from sklearn.cluster import KMeans
import warnings
from tqdm import tqdm


def get_point_reducer(point_reducer_version='cy'):
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
	return RPPointReducer


def cal_inertia(X, y_pred, centers, weight):
	n_clusters = centers.shape[0]
	inertia = 0
	for i in range(n_clusters):
		idx = np.where(y_pred == i)
		inertia += (np.square(X[idx] - centers[i]).sum(axis=1) * weight[idx]).sum()
	return inertia


def cal_cluster_variance(X, y_pred, centers, weight):
	denom = (X.shape[0] - centers.shape[0]) * X.shape[1]
	return cal_inertia(X, y_pred, centers, weight) / denom


def cal_log_likelihood(X, y_pred, centers, weight, eps=1e-100):
	log_likelihood = 0
	n_clusters = centers.shape[0]
	variance = max(eps, cal_cluster_variance(X, y_pred, centers, weight))
	total_weight = weight.sum()
	for i in range(n_clusters):
		group_size = weight[y_pred == i].sum()
		log_likelihood += group_size * np.log(group_size)
		log_likelihood -= group_size * np.log(total_weight)
		log_likelihood -= 0.5 * group_size * X.shape[1] * np.log(2.0 * np.pi * variance)
		log_likelihood -= 0.5 * X.shape[1] * (group_size - 1)
	return log_likelihood


def cal_bic(X, y_pred, centers, weight=None, weight_norm=False):
	"""Bayesian information criterion;
		Ref: Pelleg, Dan et.al., "X-means: Extending k-means with efficient estimation of the number of clusters.", 2000.
		Ref: 《Notes on Bayesian Information Criterion Calculation for X-Means Clustering》
		Ref: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans
		Ref: https://en.wikipedia.org/wiki/Bayesian_information_criterion#cite_note-Priestley-5
	Args:
		X (np.ndarray): (n_samples, n_features)
		y_pred (np.ndarray): (n_samples,)
		centers (np.ndarray): (k, n_features)
	Returns:
		float: BIC score (the higher the better)
	"""
	if weight is None:
		weight = np.ones(X.shape[0], dtype=X.dtype)
	if weight_norm:
		weight = weight * (weight.shape[0] / weight.sum())
	para_num = centers.shape[0] * (X.shape[1] + 1)
	return cal_log_likelihood(X, y_pred, centers, weight) - 0.5 * para_num * np.log(X.shape[0])


def cal_aic(X, y_pred, centers, weight=None):
	"""Ref: https://en.wikipedia.org/wiki/Akaike_information_criterion
	"""
	if weight is None:
		weight = np.ones(X.shape[0], dtype=X.dtype)
	para_num = centers.shape[0] * (X.shape[1] + 1)
	return cal_log_likelihood(X, y_pred, centers, weight) - para_num


def get_bic_optimal_k_wrapper(args):
	"""
	Args:
		args (list):
			X (np.ndarray or scipy.sparse.csr_matrix): (n_samples, n_features)
			k_range (list)
			k_repeat (int): Times to run kmeans for each k to calculate averaged bic score
			point_reducer_initializer (RPPointReducerBase)
			point_reducer_kwargs (dict): kwargs for point_reducer_initializer
	Returns:
		list: bic scores; length == len(k_range)
	"""
	scipy.random.seed()
	X, k_range, k_repeat, point_reducer_initializer, point_reducer_kwargs = args
	pr = point_reducer_initializer(**point_reducer_kwargs)
	X, weight = pr.fit_transform(X)[:-2]    # get skeleton
	bic_list = []
	for k in k_range:
		bic_k_list = []
		for repeat_id in range(k_repeat):
			while True:
				try:
					clt = KMeans(n_clusters=k)
					y_pred = clt.fit_predict(X, sample_weight=weight)
					# inertia = clt.inertia_ * (weight.sum() / X.shape[0])
				except IndexError:
					print('Index error raised, re-run kmeans for skeleton to initialize centers. This may due to bugs in sklearn.')
				else:
					break
			bic_k_list.append(cal_bic(X, y_pred, clt.cluster_centers_, weight))
		bic_list.append(np.mean(bic_k_list))
	return bic_list


def select_k_with_bic(X, kmax, kmin=2, ske_repeat=30, k_repeat=5, kneedle_s=3.0,
		point_reducer_version='cy', point_reducer_kwargs=None, n_jobs=-1):
	"""The bic score will be calculated for each k in range(kmin, kmax+1). And the optimal k will be the knee point of k-bic curve.
	Args:
		X (np.ndarray or scipy.sparse.csr_matrix):
		kmax (int)
		kmin (int): default: 2
		ske_repeat (int): default: 40
			Running times to generate skeleton
		k_repeat (int): default: 5
			Running times to calculate bic for each k in range(kmin, kmax+1) for every ske_repeat.
		kneedle_s (float): default: 5.0
			Sensitivity of Kneedle algorithm, S.
			See Satopaa, Ville, et al., "Finding a "kneedle" in a haystack: Detecting knee points in system behavior", 2011.
		point_reducer_version (str): {'cy', 'py'}, default: 'cy'
			Version of point reducer module. 'cy' is for cython version and 'py' is for python version.
			If cython version is failed to import, python version will be used instead.
		point_reducer_kwargs (dict or None): default: None
			kwargs for point reducer.
		cpu_use (int): default: -1
			The number of jobs to use for the computation. This works by computing each of the ske_repeat runs in parallel.
			-1 means using all processors.
	Returns:
		int: optimal k
		list: [bic_list(1), ..., bic_list(ske_repeat)]; bic_list(i) = [bic_score(kmin), bic_score(kmin+1), ..., bic_score(kmax)];
			bic_score(k) = np.mean([bic_score for k in range(k_repeat)]).
		list: [kmin, kmin+1, ..., kmax]
	"""
	cpu_use = cpu_count() if n_jobs == -1 else n_jobs
	point_reducer_initializer = get_point_reducer(point_reducer_version)
	k_range = list(range(kmin, kmax+1))
	point_reducer_kwargs = point_reducer_kwargs or {}
	with Pool(cpu_use) as pool:
		args_list = [(X, k_range, k_repeat, point_reducer_initializer, point_reducer_kwargs) for i in range(ske_repeat)]
		bic_lists = []
		for bic_list in tqdm(pool.imap_unordered(get_bic_optimal_k_wrapper, args_list), total=len(args_list), leave=False):
			bic_lists.append(bic_list)
	k_list = []
	s_range = ([] if int(kneedle_s) == kneedle_s else [kneedle_s]) + list(range(int(kneedle_s), 0, -1))
	for bic_list in bic_lists:
		predict_k = None
		for s in s_range:
			kl = KneeLocator(k_range, bic_list, curve='concave', direction='increasing', S=s)
			if kl.knee is not None:
				predict_k = kl.knee
				break
		assert predict_k is not None
		k_list.append(predict_k)
	optimal_k = int(round(np.mean(k_list)))
	return optimal_k, bic_lists, k_range


if __name__ == '__main__':
	pass