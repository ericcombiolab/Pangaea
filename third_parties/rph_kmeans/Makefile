default: build

FILES=point_reducer_cy setup.py README.md MANIFEST.in

point_reducer_cy:
	cython rph_kmeans/_point_reducer_cy.pyx -o rph_kmeans/_point_reducer_cy.by_cython.cpp --cplus

build: ${FILES}
	python3 setup.py build_ext -i

sdist: ${FILES}
	python3 setup.py sdist

publish: clean sdist
	twine upload dist/*

clean:
	rm -rf dist build rph_kmeans.egg-info
	rm -f rph_kmeans/*.by_cython.* *.so

