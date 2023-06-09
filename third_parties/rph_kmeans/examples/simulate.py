"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import os
# To precisely count the clustering CPU time, activate the following 3 lines when running clustering.
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import matplotlib.pyplot as plt
plt.switch_backend('agg')
plt.rcParams['axes.unicode_minus'] = False
import seaborn as sns

import scipy.sparse as sp
import itertools
from time import time, process_time
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from collections import Counter
from rph_kmeans import RPHKMeans, select_k_with_bic


def simulate_2d_gaussian_imbanlance():
	dist = 5
	sample_num_list = [5000, 100, 100, 100, 100]
	mu_list = [[0, 0], [dist, dist], [-dist, -dist], [dist, -dist], [-dist, dist]]
	cov_list = [1.0, 1.0, 1.0, 1.0, 1.0]

	y_true = np.array([i for i in range(len(sample_num_list)) for _ in range(sample_num_list[i])])
	X = []
	for sample_num, mu, cov in zip(sample_num_list, mu_list, cov_list):
		cov = np.identity(2) * cov
		X.append(np.random.multivariate_normal(mu, cov, size=sample_num))
	return np.vstack(X), y_true


def dot_plot(figpath, x, y, labels, title, sizes=None):
	"""
	Args:
		sizes (dict or list or tuple or None): {label: size}
	"""
	df = pd.DataFrame({'x': x, 'y': y, 'label': labels})
	plt.figure(figsize=(15, 15))
	ax = plt.axes()
	palette = sns.color_palette("muted", len(set(labels)))
	sns.scatterplot(x='x', y='y', hue='label', data=df, ax=ax, palette=palette,
		size='label' if sizes is not None else None, sizes=sizes)
	ax.set_title(title)
	plt.savefig(figpath)
	plt.close()


def reset_y_pred(y_pred, y_true):
	y_pred_unique = np.unique(y_pred)
	old2new = { }
	for old_lb in y_pred_unique:
		new_lb = Counter(y_true[np.where(y_pred == old_lb)]).most_common()[0][0]
		old2new[old_lb] = new_lb
	assert len(old2new.keys()) == len(y_pred_unique) and len(old2new.values()) == len(y_pred_unique)
	return np.array([old2new[old_lb] for old_lb in y_pred], dtype=np.int32)


def run_cluster(clt_initializer, clt_kwargs, X, y_true, figpath, clt_name, repeat=10, draw=True):
	p_time, t_time = process_time(), time()
	result_dicts = []
	for repeat_id in range(repeat):
		clt = clt_initializer(**clt_kwargs)
		y_pred = clt.fit_predict(X)
		p_time_spend, t_time_spend = process_time() - p_time, time() - t_time
		ari, nmi = adjusted_rand_score(y_true, y_pred), normalized_mutual_info_score(y_true, y_pred),
		result_dicts.append({
			'CLUSTER_REAL_TIME': t_time_spend,
			'CLUSTER_CPU_TIME': p_time_spend,
			'ARI': ari,
			'NMI': nmi
		})
	print('{}: \n{}'.format(clt_name, combine_metric_dicts(result_dicts)))
	if draw:
		dot_plot(figpath, X[:, 0], X[:, 1], y_pred, title='{} (ARI={:.2f}; NMI={:.2f})'.format(clt_name, ari, nmi))


def combine_metric_dicts(dlist):
	"""
	Args:
		dlist (list): [{metric_name: score}, ...]
	Returns:
		dict: {metric_name: score_mean (score_std)}
	"""
	d = {}
	for k in dlist[0]:
		score_list = [metric_dict[k] for metric_dict in dlist]
		ave = np.mean(score_list)
		std = np.std(score_list)
		d[k] = '{:.3f} ({:.3f})'.format(ave, std)
	return d


def test_performance(dense_X, y_true):
	fig_save_folder = 'performance_test'
	os.makedirs(fig_save_folder, exist_ok=True)
	sparse_X = sp.csr_matrix(dense_X)
	dot_plot(os.path.join(fig_save_folder, 'y_true'), dense_X[:, 0], dense_X[:, 1], y_true, 'y_true')

	n_clusters = len(np.unique(y_true))
	point_reducer_versions = ['py', 'cy']
	is_sparse_list = [False, True]
	bkt_improve_dicts = [
		{'bkt_improve': None},
		{'bkt_improve': 'radius', 'radius_divide':0.5},
		{'bkt_improve': 'min_bkt_size', 'bkt_size_keepr':0.8},
		{'bkt_improve': 'min_center_dist', 'center_dist_keepr':0.8},
	]

	clt_initializer = RPHKMeans
	for prv, is_sparse, bkt_imp_dict in itertools.product(point_reducer_versions, is_sparse_list, bkt_improve_dicts):
		X = sparse_X if is_sparse else dense_X
		clt_name = 'rph_kmeans_{}_{}_imp_{}'.format(prv, 'sparsex' if is_sparse else 'densex',
			bkt_imp_dict['bkt_improve'])
		print('\nrunning {} -------------------------'.format(clt_name))
		clt_kwargs = {**{'n_clusters':n_clusters, 'point_reducer_version':prv, 'verbose':True}, **bkt_imp_dict}
		run_cluster(clt_initializer, clt_kwargs, X, y_true, os.path.join(fig_save_folder, clt_name), clt_name, draw=False)

	clt_initializer = KMeans
	for init, n_init in itertools.product(['k-means++', 'random'], [1, 5, 10]):
		clt_name = 'kmeans_{}_{}'.format(init, n_init)
		print('\nrunning {} -------------------------'.format(clt_name))
		clt_kwargs = {'n_clusters':n_clusters, 'init':init, 'n_init':n_init}
		run_cluster(clt_initializer, clt_kwargs, dense_X, y_true, os.path.join(fig_save_folder, clt_name), clt_name, draw=False)


def show_pipeline(dense_X, y_true):
	fig_save_folder = 'pipeline_draw'
	os.makedirs(fig_save_folder, exist_ok=True)
	n_clusters = len(np.unique(y_true))
	dot_plot(os.path.join(fig_save_folder, 'y_true.png'), dense_X[:, 0], dense_X[:, 1], y_true, 'True label')

	# RPH-KMeans =======================================================================
	clt = RPHKMeans(n_clusters=n_clusters)
	y_pred = clt.fit_predict(dense_X)
	y_pred = reset_y_pred(y_pred, y_true)   # To make the picture consistent

	draw_points = np.vstack([clt.reduced_X_, clt.init_centers_])
	draw_labels = ['reduced point'] * clt.reduced_X_.shape[0] + ['initial center'] * n_clusters
	dot_plot(fig_save_folder + os.sep + 'rph_kmeans_reduced_points.png', draw_points[:, 0], draw_points[:, 1],
		draw_labels, 'Reduced points and initial centers', sizes={'reduced point': 40, 'initial center': 100})

	draw_points = np.vstack([dense_X, clt.cluster_centers_])
	draw_labels = ['data point'] * dense_X.shape[0] + ['cluster center'] * n_clusters
	dot_plot(fig_save_folder + os.sep + 'rph_kmeans_cluster_centers.png', draw_points[:, 0], draw_points[:, 1],
		draw_labels, 'Data points and cluster centers', sizes={'data point': 40, 'cluster center': 100})

	ari, nmi = adjusted_rand_score(y_true, y_pred), normalized_mutual_info_score(y_true, y_pred),
	dot_plot(fig_save_folder + os.sep + 'rph_kmeans_y_pred.png',
		dense_X[:, 0], dense_X[:, 1], y_pred, title='rph-kmeans (ARI={:.2f}; NMI={:.2f})'.format(ari, nmi))

	# KMeans (kmeans++ init) =======================================================================
	clt = KMeans(n_clusters=n_clusters, n_init=1, init='k-means++')
	y_pred = clt.fit_predict(dense_X)
	# y_pred = reset_y_pred(y_pred, y_true)  # To make the picture consistent

	draw_points = np.vstack([dense_X, clt.cluster_centers_])
	draw_labels = ['data point'] * dense_X.shape[0] + ['cluster center'] * n_clusters
	dot_plot(fig_save_folder + os.sep + 'kmeans(kmeans++)_cluster_centers.png', draw_points[:, 0], draw_points[:, 1],
		draw_labels, 'Data points and cluster centers', sizes={'data point':40, 'cluster center':100})

	ari, nmi = adjusted_rand_score(y_true, y_pred), normalized_mutual_info_score(y_true, y_pred),
	dot_plot(fig_save_folder + os.sep + 'kmeans(kmeans++)_y_pred.png',
		dense_X[:, 0], dense_X[:, 1], y_pred, title='kmeans(kmeans++) (ARI={:.2f}; NMI={:.2f})'.format(ari, nmi))

	# KMeans (random init) =======================================================================
	clt = KMeans(n_clusters=n_clusters, n_init=1, init='random')
	y_pred = clt.fit_predict(dense_X)
	# y_pred = reset_y_pred(y_pred, y_true)  # To make the picture consistent

	draw_points = np.vstack([dense_X, clt.cluster_centers_])
	draw_labels = ['data point'] * dense_X.shape[0] + ['cluster center'] * n_clusters
	dot_plot(fig_save_folder + os.sep + 'kmeans(random)_cluster_centers.png', draw_points[:, 0], draw_points[:, 1],
		draw_labels, 'Data points and cluster centers', sizes={'data point':40, 'cluster center':100})

	ari, nmi = adjusted_rand_score(y_true, y_pred), normalized_mutual_info_score(y_true, y_pred),
	dot_plot(fig_save_folder + os.sep + 'kmeans(random)_y_pred.png',
		dense_X[:, 0], dense_X[:, 1], y_pred, title='kmeans(random) (ARI={:.2f}; NMI={:.2f})'.format(ari, nmi))


def test_k_selection(dense_X, y_true):
	fig_save_folder = 'k_selection'
	os.makedirs(fig_save_folder, exist_ok=True)
	n_clusters_true = len(np.unique(y_true))

	kmax = n_clusters_true * 3
	n_clusters_pred, bic_lists, k_range = select_k_with_bic(dense_X, kmax=kmax, point_reducer_kwargs={'verbose': 0})
	print(f'True cluster number = {n_clusters_true}; Predicted cluster number = {n_clusters_pred}')

	plt.figure(); ax = plt.axes()
	bic_list_draw = bic_lists[0]
	sns.lineplot(x='K', y='BIC', data=pd.DataFrame({'K': k_range, 'BIC': bic_list_draw}), ax=ax)
	v_min, v_max = min(bic_list_draw), max(bic_list_draw)
	plt.vlines(n_clusters_true, v_min, v_max, colors='k', label=f'True k', linestyles='--')
	plt.vlines(n_clusters_pred, v_min, v_max, colors='r', label=f'Predicted k', linestyles='--')
	plt.legend()
	plt.savefig(os.path.join(fig_save_folder, 'K-BIC.png'))
	plt.close()


if __name__ == '__main__':
	dense_X, y_true = simulate_2d_gaussian_imbanlance()

	test_performance(dense_X, y_true)
	show_pipeline(dense_X, y_true)
	test_k_selection(dense_X, y_true)


