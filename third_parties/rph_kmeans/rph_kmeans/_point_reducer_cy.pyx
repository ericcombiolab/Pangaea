# cython: profile=False
# cython: boundscheck=False, wraparound=False, cdivision=True

"""
@Author: Yu Huang
@Email: yuhuang-cst@foxmail.com
"""

import numpy as np
cimport numpy as np
cimport cython
from cython cimport floating

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INT
ctypedef np.uint32_t UINT
ctypedef unsigned int uint

np.import_array()

# python3 setup.py build_ext --inplace

cdef extern from "_point_reducer_cy_lib.h":
    uint get_ary_labels_(int* ary, uint* labels, uint r, uint c)
    void update_labels_(uint* labels, uint* bkt_ary, uint labels_r)
    void update_densex_and_weight_[T](T* X, T* weight, uint* bkt_ary, uint r, uint c, uint bkt_num, T* out_X, T* out_weight)
    void update_sparsex_and_weight_[T](T* data, int* indices, int* indptr, T* weight, uint* bkt_ary, uint r, uint c, uint bkt_num, T* out_X, T* out_weight)
    uint densex_radius_bkt_improve_[T](T* X, uint* bkt_ary, T R2, uint r, uint c, uint bkt_num)
    uint sparsex_radius_bkt_improve_[T](T* data, int* indices, int* indptr, uint* bkt_ary, T R2, uint r, uint c, uint bkt_num)


def get_ary_labels(np.ndarray[int, ndim=2, mode='c'] ary):
    cdef np.ndarray[uint, ndim=1] labels
    cdef np.ndarray[uint, ndim=1, mode='c'] bkt_ary
    bkt_ary = np.zeros((ary.shape[0],), dtype=np.uint32)
    bkt_num = get_ary_labels_(&ary[0,0], &bkt_ary[0], ary.shape[0], ary.shape[1])
    return bkt_ary, bkt_num


def update_labels(np.ndarray[uint, ndim=1, mode='c'] labels, np.ndarray[uint, ndim=1, mode='c'] bkt_ary):
    update_labels_(&labels[0], &bkt_ary[0], labels.shape[0])




def update_densex_and_weight(
        np.ndarray[floating, ndim=2, mode='c'] X, np.ndarray[floating, ndim=1, mode='c'] weight,
		np.ndarray[uint, ndim=1, mode='c'] bkt_ary, uint bkt_num):
    cdef uint row_num = X.shape[0]
    cdef uint col_num = X.shape[1]
    cdef np.ndarray[floating, ndim=2] out_X
    cdef np.ndarray[floating, ndim=1] out_weight
    dtype = np.float32 if floating is float else np.float64
    out_X = np.zeros((bkt_num, col_num), dtype=dtype)
    out_weight = np.zeros((bkt_num,), dtype=dtype)

    update_densex_and_weight_(&X[0, 0], &weight[0], &bkt_ary[0], row_num, col_num, bkt_num, &out_X[0, 0], &out_weight[0])
    return out_X, out_weight


def update_sparsex_and_weight(X, np.ndarray[floating, ndim=1, mode='c'] weight, np.ndarray[uint, ndim=1, mode='c'] bkt_ary, uint bkt_num):
    cdef uint row_num = X.shape[0]
    cdef uint col_num = X.shape[1]
    cdef np.ndarray[floating, ndim=1] data = X.data
    cdef np.ndarray[int, ndim=1] indices = X.indices
    cdef np.ndarray[int, ndim=1] indptr = X.indptr
    cdef np.ndarray[floating, ndim=2] out_X
    cdef np.ndarray[floating, ndim=1] out_weight
    dtype = np.float32 if floating is float else np.float64
    out_X = np.zeros((bkt_num, col_num), dtype=dtype)
    out_weight = np.zeros((bkt_num,), dtype=dtype)

    update_sparsex_and_weight_(&data[0], &indices[0], &indptr[0], &weight[0], &bkt_ary[0], row_num, col_num, bkt_num, &out_X[0, 0], &out_weight[0])
    return out_X, out_weight


def densex_radius_bkt_improve(np.ndarray[floating, ndim=2, mode='c'] X, np.ndarray[uint, ndim=1, mode='c'] bkt_ary, floating R2, uint bkt_num):
    cdef uint row_num = X.shape[0]
    cdef uint col_num = X.shape[1]
    return densex_radius_bkt_improve_(&X[0, 0], &bkt_ary[0], R2, row_num, col_num, bkt_num)


def sparsex_radius_bkt_improve(
        np.ndarray[floating, ndim=1, mode='c'] data, np.ndarray[int, ndim=1, mode='c'] indices, np.ndarray[int, ndim=1, mode='c'] indptr,
        uint row_num, uint col_num, np.ndarray[uint, ndim=1, mode='c'] bkt_ary, floating R2, uint bkt_num):
    return sparsex_radius_bkt_improve_(&data[0], &indices[0], &indptr[0], &bkt_ary[0], R2, row_num, col_num, bkt_num)
