set -e
# install jellyfish
cd third_parties/; tar -xf jellyfish-2.3.0.tar.gz; cd -
cd third_parties/jellyfish-2.3.0/
./configure --prefix=$PWD/../jellyfish_pkg
make -j 4   
make install
cd -

# install conda env pangaea
conda env create -f environment.yaml
source $CONDA_PREFIX/bin/activate pangaea
if [ $CONDA_DEFAULT_ENV == "pangaea" ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea activated"
else
    echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea not activated"
    exit
fi
# install rph_kmeans
cd third_parties/rph_kmeans
python3 setup.py install
cd -

# build util binaries
cd cpptools && make && cd -