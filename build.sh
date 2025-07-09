#!/bin/bash
set -e

# install conda env pangaea
mamba env create -f environment.yaml && conda activate pangaea

# install athena-meta
mamba create -n athena-meta bioconda::athena_meta

# install pytorch according to https://pytorch.org/get-started/locally/
# for example, if cpu only
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# install rph_kmeans
conda activate pangaea
cd third_parties/rph_kmeans
python setup.py install
cd ../../

# build cpp bins
cd src/cpptools
if [ -d "build" ]; then
  rm -rf build
fi
mkdir build && cd build
cmake ..
make
cd ../../../