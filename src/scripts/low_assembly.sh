#!/bin/bash
assembly_dir=$1
contigs=$2
covcut=$3
thread=$4
assembler=$5
BINDIR=$(dirname "$0")

set -e
if [ $assembler == "spades" ];then
    if [ ! -d $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades ]; then
        metaspades.py --12 $assembly_dir/contigs.megahit_cut$covcut.low_abd.fq --only-assembler -m 5000 -t $thread --untrusted-contigs $contigs -o $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades > /dev/null
    fi
    if [ ! -f $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades/contigs.fasta ]; then
        metaspades.py --continue -o $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades > /dev/null
    fi
else
    if [ ! -f $assembly_dir/contigs.megahit_cut$covcut.low_abd.megahit/final.contigs.fa ]; then
        megahit --12 $assembly_dir/contigs.megahit_cut$covcut.low_abd.fq -t $thread -o $assembly_dir/contigs.megahit_cut$covcut.low_abd.megahit > /dev/null 2>&1
    fi
fi