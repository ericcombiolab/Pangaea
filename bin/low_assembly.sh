#!/bin/bash
assembly_dir=$1
contigs=$2
covcut=$3
thread=$4
BINDIR=$(dirname "$0")

set -e

if [ ! -d $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades ]; then
    metaspades.py --12 $assembly_dir/contigs.megahit_cut$covcut.low_abd.fq --only-assembler -m 5000 -t $thread --untrusted-contigs $contigs -o $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades > /dev/null
fi
if [ ! -f $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades/contigs.fasta ]; then
    metaspades.py --continue -o $assembly_dir/contigs.megahit_cut$covcut.low_abd.spades > /dev/null
fi