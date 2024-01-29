#!/bin/bash
# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.tar
# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vOct22_CHOCOPhlAnSGB_202212_bt2.md5
# wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar
threads=`nproc`
input_reads=$1
# this file's path
script_path=$(dirname "$0")
mkdir -p ./metaphlan_tmp
metaphlan $input_reads --input_type fastq --bowtie2db $DATA/metaphlan/ --nproc $threads --bowtie2out metaphlan_tmp/metagenome_from_reads.bowtie2.bz2 -o metaphlan_tmp/profiled.txt
python $script_path/metaphlan_tables.py  metaphlan_tmp/profiled.txt metaphlan_tmp/profiled.txt > profiles_table.tsv
Rscript $script_path/calculate_diversity.R -d alpha -m shannon -f  profiles_table.tsv -o metaphlan_tmp/diversity_analysis
