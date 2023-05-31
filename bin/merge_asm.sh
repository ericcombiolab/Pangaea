#!/bin/bash
cluster_dir=$1
assembly_dir=$2
local_assembly=$3
athena=$4
contigs=$5
snakefile=$6
BINDIR=$(dirname "$0")

set -e

cat $assembly_dir/*.spades/contigs.fasta $cluster_dir/contigs.megahit.fa $local_assembly > $assembly_dir/contigs.low_abd.binning.local.fa
$BINDIR/parse_header $assembly_dir/contigs.low_abd.binning.local.fa contig_ > $assembly_dir/contigs.low_abd.binning.local.renamed.fa
mv $assembly_dir/contigs.low_abd.binning.local.renamed.fa $assembly_dir/contigs.low_abd.binning.local.fa

$BINDIR/merge_olc.py $contigs $assembly_dir/contigs.low_abd.binning.local.fa $assembly_dir/contigs.low_abd.binning.local.asm

if [ ! -d $assembly_dir/quickmerge ]; then
    mkdir $assembly_dir/quickmerge
fi

athena=`realpath $athena`
cd $assembly_dir/quickmerge

# quickmerge
merge_wrapper.py ../contigs.low_abd.binning.local.asm/final.asm.fa $athena
$BINDIR/parse_header merged_out.fasta contig_ > merged_out.renamed.fasta
mv merged_out.renamed.fasta merged_out.fasta

echo 'sample_name: "circular"' > circular.yaml
echo 'contigs: "merged_out.fasta"' >> circular.yaml
echo 'reads: "../contigs.low_abd.binning.local.fa"' >> circular.yaml
{
    snakemake -s $snakefile --configfile circular.yaml --cores 100 > /dev/null 2>&1 &&
    cp circular/3.circularization/4.circular_circularized.fasta ../../final.asm.fa
} || {
    cp merged_out.fasta ../../final.asm.fa
}
