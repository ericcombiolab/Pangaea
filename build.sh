# install libs
cd third_parties/jelleyfish-2.3.0
./configure --prefix=$PWD/../jellyfish_pkg
make -j 4   
make install
cd -

# install conda env 
conda env create -f environment.yaml
conda activate pangaea

# install rph_kmeans
cd third_parties/rph_kmeans
python3 setup.py install
cd -

# build util binaries
cd Pangaea/cpptools && make && cd -