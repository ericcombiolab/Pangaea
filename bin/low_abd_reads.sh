#!/bin/bash
cluster_dir=$1
assembly_dir=$2
covcut=$3
BINDIR=$(dirname "$0")

set -e

if [ ! -f $assembly_dir/contigs.megahit_cut$covcut.low_abd.fq ]; then
    $BINDIR/extract_unmapped -b $assembly_dir/contigs.megahit.name_sorted.bam -c $assembly_dir/contigs.megahit.depth -f $covcut -o $assembly_dir/contigs.megahit_cut$covcut
    seqtk subseq $cluster_dir/contigs.megahit.fa $assembly_dir/contigs.megahit_cut$covcut.list > $assembly_dir/contigs.megahit_cut$covcut.high_abd.fa
fi