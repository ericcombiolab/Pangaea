#!/bin/bash
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

conda install -c bioconda metaphlan

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