CXX = g++
CPPFLAGS = -I. -O3
LDFLAGS = -L. -lgzstream -lz -lpthread
AR = ar cr

default: count_kmer count_tnf extract_reads extract_unmapped parse_header
	cp count_kmer count_tnf extract_reads extract_unmapped parse_header ../bin

gzstream.o : gzstream.C gzstream.h
	${CXX} ${CPPFLAGS} -c -o gzstream.o gzstream.C
libgzstream.a : gzstream.o
	${AR} libgzstream.a gzstream.o
count_kmer : count_kmer.cpp cmdline.h libgzstream.a
	${CXX} ${CPPFLAGS} -o count_kmer count_kmer.cpp ${LDFLAGS}
count_tnf : count_tnf.cpp cmdline.h libgzstream.a
	${CXX} ${CPPFLAGS} -o count_tnf count_tnf.cpp ${LDFLAGS}
extract_reads : extract_reads.cpp cmdline.h libgzstream.a
	${CXX} ${CPPFLAGS} -o extract_reads extract_reads.cpp ${LDFLAGS}
down_reads : down_reads.cpp cmdline.h libgzstream.a
	${CXX} ${CPPFLAGS} -o down_reads down_reads.cpp ${LDFLAGS}
extract_unmapped : extract_unmapped.cpp cmdline.h
	${CXX} ${CPPFLAGS} -o extract_unmapped extract_unmapped.cpp -lhts -fopenmp
parse_header: parse_header.cpp
	${CXX} ${CPPFLAGS} -o parse_header parse_header.cpp

clean:
	rm -f *.o libgzstream.a count_kmer count_tnf extract_reads extract_unmapped
	rm -f ../bin/count_kmer ../bin/count_tnf ../bin/extract_reads ../bin/extract_unmapped