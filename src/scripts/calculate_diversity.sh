#!/bin/bash
# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.md5
# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar
threads=`nproc`
input_reads=$1
if [ ! -f $input_reads ];then
    echo "Error: $input_reads not found"
    exit 1
fi

find_latest_mpa_index() {
    local database_dir="$1"

    if [[ ! -d "$database_dir" ]]; then
        echo "Error: The directory '$database_dir' does not exist." >&2
        return 1
    fi

    # find all MetaPhlAn database index files
    local index_files=($(ls "$database_dir" | grep -E '^(mpa_v[^ ]+)\.1\.bt2l$' | grep -v '\.rev\.'))
    if [[ ${#index_files[@]} -eq 0 ]]; then
        echo "Error: No valid MetaPhlAn database index found in '$database_dir'" >&2
        return 1
    fi

    # extract the version prefixes from the index files
    local versions=()
    for file in "${index_files[@]}"; do
        # eg. mpa_vOct22_CHOCOPhlAnSGB_202212
        local prefix=$(echo "$file" | sed -E 's/\.1\.bt2l$//')
        versions+=("$prefix")
    done

    local latest_version=$(printf '%s\n' "${versions[@]}" | sort | tail -n 1)

    echo "$latest_version"
}


# this file's path
script_path=$(dirname "$0")
metaphlan_db=$2
cluster_path=$3
# check the validity of all input parameters
if [ -z "$input_reads" ] || [ -z "$metaphlan_db" ] || [ -z "$cluster_path" ]; then
    echo "Usage: $0 <input_reads> <metaphlan_db> <cluster_path>"
    exit 1
fi

if [ ! -d $cluster_path ] || [ ! -d $metaphlan_db ] || [ ! -d $script_path ]; then
    echo "Error: Invalid directory path(s) provided."
    exit 1
fi

mkdir -p $cluster_path/metaphlan_tmp
latest_index=$(find_latest_mpa_index $metaphlan_db)

metaphlan $input_reads --offline --index $latest_index --input_type fastq --bowtie2db $metaphlan_db --nproc $threads --bowtie2out $cluster_path/metaphlan_tmp/metagenome_from_reads.bowtie2.bz2 -o $cluster_path/metaphlan_tmp/profiled.txt

python $script_path/metaphlan_tables.py  $cluster_path/metaphlan_tmp/profiled.txt $cluster_path/metaphlan_tmp/profiled.txt > $cluster_path/metaphlan_tmp/profiles_table.tsv
mkdir -p $cluster_path/metaphlan_tmp/diversity_analysis

Rscript $script_path/calculate_diversity.R -d alpha -m shannon -f  $cluster_path/metaphlan_tmp/profiles_table.tsv -o $cluster_path/metaphlan_tmp/diversity_analysis
