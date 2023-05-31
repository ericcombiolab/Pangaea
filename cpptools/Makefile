LDFLAGS = -L./lib -lpthread -lgzstream -lz -Wl,-rpath $(PWD)/lib

default: count_kmer count_tnf extract_reads extract_unmapped parse_header

count_kmer : count_kmer.cpp util.cpp
	g++ -O3 -o ../bin/count_kmer count_kmer.cpp util.cpp -I./include  ${LDFLAGS}
count_tnf : count_tnf.cpp
	g++ -O3 -o ../bin/count_tnf count_tnf.cpp -I./include ${LDFLAGS}
extract_reads : extract_reads.cpp util.cpp
	g++ -O3 -o ../bin/extract_reads extract_reads.cpp util.cpp -I./include ${LDFLAGS}
extract_reads_long : extract_reads_longreads.cpp util.cpp
	g++ -O3 -o ../bin/extract_reads_long extract_reads_longreads.cpp util.cpp -I./include ${LDFLAGS}
extract_unmapped : extract_unmapped.cpp
	g++ -O3 -o ../bin/extract_unmapped extract_unmapped.cpp -fopenmp -I./include -lhts -lbz2 ${LDFLAGS}
extract_unmapped_long : extract_unmapped_longreads.cpp
	g++ -O3 -o ../bin/extract_unmapped_longreads extract_unmapped_longreads.cpp -fopenmp -I./include -lhts -lbz2 ${LDFLAGS}
parse_header: parse_header.cpp
	g++ -O3 -o ../bin/parse_header parse_header.cpp -I./include ${LDFLAGS}
add_barcode : add_barcode.cpp
	g++ -O3 -o ../bin/add_barcode add_barcode.cpp -I./include -lhts -fopenmp
assign_barcodes: assign_barcodes.cpp
	g++ -O3 -o ../bin/assign_barcodes assign_barcodes.cpp -I./include -fopenmp ${LDFLAGS}

clean:
	rm -f ../bin/count_kmer ../bin/count_tnf ../bin/extract_reads ../bin/extract_unmapped