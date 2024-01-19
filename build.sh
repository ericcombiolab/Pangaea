set -e

while getopts "c:" opt; do
    case $opt in
        d)
            metapglan4_DB=$OPTARG
            ;;
        ?)
            echo "Usage: build.sh [-d metapglan4_DB]"
            exit 1
            ;;
    esac
done


# install jellyfish
cd third_parties/; tar -xf jellyfish-2.3.0.tar.gz; cd -
cd third_parties/jellyfish-2.3.0/
./configure --prefix=$PWD/../jellyfish_pkg
make -j 4   
make install
cd -

# install rph_kmeans
cd third_parties/rph_kmeans
python3 setup.py install
cd -

# build cpp bins
cd cpptools && make && cd -

# install conda env pangaea
conda env create -f environment.yaml

source $CONDA_PREFIX/bin/activate pangaea
if [ $CONDA_DEFAULT_ENV == "pangaea" ];then
    echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea activated"
else
    echo "`date "+%Y-%m-%d %H:%M:%S"` pangaea not activated"
    exit
fi

# if metaphlan4_DB is not specified by user, use the default one
if [ -z $metapglan4_DB ];then
    metapglan4_DB=$PWD/metaphlan4_DB
fi
fi

# install metaphlan4 Database
if [ ! -d "$metapglan4_DB" ];then
    mkdir -p $metapglan4_DB
fi
metaphlan --install --bowtie2db $metapglan4_DB
echo "`date "+%Y-%m-%d %H:%M:%S"` metaphlan4_DB installed into $metapglan4_DB"