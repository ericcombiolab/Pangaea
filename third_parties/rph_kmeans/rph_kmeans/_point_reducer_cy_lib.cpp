#include <string>
#include <unordered_map>
#include "_point_reducer_cy_lib.h"

uint get_ary_labels_(int* ary, uint* bkt_ary, uint r, uint c){
	int row_size = sizeof(int) * c;
	unordered_map< string, list<int> > m;
	char *p;
	string row_bytes(row_size, '\0');
	for (uint i = 0; i < r; ++i){
		p = (char*) &ary[i*c];
		row_bytes.assign(p, row_size);
		if (m.find(row_bytes) == m.end()){ // not found
			list<int> ll = list<int>();
			ll.push_back(i);
			m[row_bytes] = ll;
		} else {
			m[row_bytes].push_back(i);
		}
	}
	uint lb = 0;
	for (auto it = m.begin(); it != m.end(); ++it){
		for (auto row_id_it = it->second.begin(); row_id_it != it->second.end(); ++row_id_it)
			bkt_ary[*row_id_it] = lb;
		lb += 1;
	}
	return lb;
}


void update_labels_(uint* labels, uint* bkt_ary, uint labels_r){
    for (uint i = 0; i < labels_r; ++i)
        labels[i] = bkt_ary[labels[i]];
}
