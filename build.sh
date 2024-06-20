#!/bin/bash
set -e

# install conda env pangaea
conda env create -f environment.yaml

# install spades(metaspades)
mamba install spades metaphlan

# install rph_kmeans
cd third_parties/rph_kmeans
python3 setup.py install
cd -

# build cpp bins
cd cpptools && make && cd -


source $CONDA_PREFIX/bin/activate pangaea
if [ $CONDA_DEFAULT_ENV == "pangaea" ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea activated"
else
    echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea not activated"
    exit
fi
