"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import scipy.sparse as sp
import numpy as np

def check_return(attrCollector):
	def outerWrapper(func):
		def wrapper(cls, *args, **kwargs):
			coll = getattr(cls, attrCollector, None)
			if coll is not None:
				return coll
			coll = func(cls, *args, **kwargs)
			setattr(cls, attrCollector, coll)
			return coll
		return wrapper
	return outerWrapper


def cal_dist2_ary_sparse(X, Y):
	if X.shape[0] != Y.shape[0]:
		Y = sp.vstack([Y] * X.shape[0])
	return (X - Y).power(2).sum(axis=1).A1


def cal_dist2_ary_dense(X, Y):
	return np.square(X - Y).sum(axis=1)


def cal_weighted_ave_sparse_vec(X, weight):
	return X.multiply(weight).sum(axis=0).A1


def cal_weighted_ave_dense_vec(X, weight):
	return np.multiply(X, weight).sum(axis=0)



