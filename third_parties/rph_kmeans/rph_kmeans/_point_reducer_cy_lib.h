#ifndef RP_KMEANS_UNIT_H
#define RP_KMEANS_UNIT_H

#include <list>
#include <limits>
using namespace std;

typedef unsigned int uint;

uint get_ary_labels_(int* ary, uint* labels, uint r, uint c);
void update_labels_(uint* labels, uint* bkt_ary, uint labels_r);
void test_get_ary_labels_();
void test_update_densex_and_weight();

template <typename T>
inline T dense_euclidian_dist2(T* ary1, T* ary2, uint size){
    T dst = 0, tmp = 0;
    for (uint i = 0; i < size; ++i){
        tmp = ary1[i] - ary2[i];
        dst += tmp * tmp;
    }
    return dst;
}

template <typename T>
inline T sparse_euclidian_dist2(T* ary1, int* indices1, int size1, T* ary2, int* indices2, int size2){
    int i = 0, j = 0;
    T dst = 0, tmp;
    while (true){
        if (i == size1){
            for (; j < size2; ++j)
                dst += ary2[j] * ary2[j];
            break;
        }
        if (j == size2){
            for (; i < size1; ++i)
                dst += ary1[i] * ary1[i];
            break;
        }
        if (indices1[i] < indices2[j]){
            dst += ary1[i] * ary1[i]; ++i;
        } else if (indices1[i] > indices2[j]){
            dst += ary2[j] * ary2[j]; ++j;
        } else {
            tmp = ary1[i] - ary2[j];
            dst += tmp * tmp; ++i; ++j;
        }
    }
    return dst;
}

template <typename T>
void update_densex_and_weight_(T* X, T* weight, uint* bkt_ary, uint r, uint c, uint bkt_num, T* out_X, T* out_weight){
	uint i, j, b, bkt_id, out_b;
	for (i = 0; i < r; ++i){
		bkt_id = bkt_ary[i];
		out_weight[bkt_id] += weight[i];
		b = i * c;
		out_b = bkt_id * c;
		for (j = 0; j < c; ++j)
			out_X[out_b + j] += X[b + j] * weight[i];
	}
	for (i = 0; i < bkt_num; ++i){
		b = i * c;
		for (j = 0; j < c; ++j)
			out_X[b + j] /= out_weight[i];
	}
}

template <typename T>
void update_sparsex_and_weight_(T* data, int* indices, int* indptr, T* weight, uint* bkt_ary, uint r, uint c, uint bkt_num, T* out_X, T* out_weight){
    uint i, j, bkt_id, b, out_b;
    int ind;
    for (i = 0; i < r; ++i){
        bkt_id = bkt_ary[i];
        out_weight[bkt_id] += weight[i];
        out_b = bkt_id * c;
        for (ind = indptr[i]; ind < indptr[i+1]; ++ind){
            j = indices[ind];
            out_X[out_b + j] += data[ind] * weight[i];
        }
    }
    for (i = 0; i < bkt_num; ++i){
		b = i * c;
		for (j = 0; j < c; ++j)
			out_X[b + j] /= out_weight[i];
	}
}

template <typename T>
uint densex_radius_bkt_improve_(T* X, uint* bkt_ary, T R2, uint r, uint c, uint bkt_num){
    list<uint>* group_centers = new list<uint>[bkt_num];
    list<uint>* g_ptr;
    T* b_ptr;
    T dst2, min_dist2;
    uint lb = 0, i, belong_to;
    for (i = 0; i < r; ++i){
        g_ptr = &group_centers[bkt_ary[i]];
        if (g_ptr->empty()){
            g_ptr->push_back(i);
            bkt_ary[i] = lb; lb += 1;
            continue;
        }
        b_ptr = X + i * c;
        min_dist2 = numeric_limits<int>::max();
        for (auto it = g_ptr->begin(); it != g_ptr->end(); ++it){
            dst2 = dense_euclidian_dist2(b_ptr, X + (*it) * c, c);
            if (dst2 < min_dist2){
                min_dist2 = dst2;
                belong_to = *it;
            }
        }
        if (min_dist2 < R2)
            bkt_ary[i] = bkt_ary[belong_to];
        else {
            g_ptr->push_back(i);
            bkt_ary[i] = lb; lb += 1;
        }
    }
    delete []group_centers;
    return lb;
}

template <typename T>
uint sparsex_radius_bkt_improve_(T* data, int* indices, int* indptr, uint* bkt_ary, T R2, uint r, uint c, uint bkt_num){
    list<uint>* group_centers = new list<uint>[bkt_num];
    list<uint>* g_ptr;
    T dst2, min_dist2;
    uint lb = 0, i, j, belong_to;
    for (i = 0; i < r; ++i){
        g_ptr = &group_centers[bkt_ary[i]];
        if (g_ptr->empty()){
            g_ptr->push_back(i);
            bkt_ary[i] = lb; lb += 1;
            continue;
        }
        min_dist2 = numeric_limits<int>::max();
        for (auto it = g_ptr->begin(); it != g_ptr->end(); ++it){
            j = *it;
            dst2 = sparse_euclidian_dist2(
                data+indptr[i], indices+indptr[i], indptr[i+1] - indptr[i],
                data+indptr[j], indices+indptr[j], indptr[j+1] - indptr[j]);
            if (dst2 < min_dist2){
                min_dist2 = dst2;
                belong_to = *it;
            }
        }
        if (min_dist2 < R2)
            bkt_ary[i] = bkt_ary[belong_to];
        else {
            g_ptr->push_back(i);
            bkt_ary[i] = lb; lb += 1;
        }
    }
    delete []group_centers;
    return lb;
}

#endif


