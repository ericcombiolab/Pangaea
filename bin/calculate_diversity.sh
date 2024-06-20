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

# this file's path
script_path=$(dirname "$0")
mkdir -p ./metaphlan_tmp
metaphlan_db=$2
if [ ! -d $metaphlan_db ];then
    echo "Not provide the metaphlan database path, use the default one"
    metaphlan $input_reads --input_type fastq --nproc $threads --bowtie2out metaphlan_tmp/metagenome_from_reads.bowtie2.bz2 -o metaphlan_tmp/profiled.txt
else
    metaphlan $input_reads --input_type fastq --bowtie2db $metaphlan_db --nproc $threads --bowtie2out metaphlan_tmp/metagenome_from_reads.bowtie2.bz2 -o metaphlan_tmp/profiled.txt
fi
python $script_path/metaphlan_tables.py  metaphlan_tmp/profiled.txt metaphlan_tmp/profiled.txt > profiles_table.tsv
Rscript $script_path/calculate_diversity.R -d alpha -m shannon -f  profiles_table.tsv -o metaphlan_tmp/diversity_analysis
