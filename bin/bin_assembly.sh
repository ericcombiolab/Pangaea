#!/bin/bash
cluster_dir=$1
assembly_dir=$2
thread=$3
reads1=$4
reads2=$5
BINDIR=$(dirname "$0")
num_jobs="\j"

set -e

if [ ! -d $assembly_dir ]; then
    mkdir $assembly_dir
fi

if [ ! -f $assembly_dir/contigs.megahit.name_sorted.bam ]; then
    if [ ! -f $cluster_dir/contigs.megahit.fa ]; then
        for fq in $cluster_dir/*.fq
        do
            while (( ${num_jobs@P} >= 5 )); do
                wait -n
            done
            if [ ! -d ${fq%.fq}.megahit ]; then
                megahit --12 $fq -t $thread -o ${fq%.fq}.megahit > /dev/null 2>&1 &
            fi
        done
        wait
        cat $cluster_dir/*.megahit/final.contigs.fa > $cluster_dir/contigs.megahit.fa
        $BINDIR/parse_header $cluster_dir/contigs.megahit.fa contig_ > $cluster_dir/contigs.megahit.renamed.fa
        mv $cluster_dir/contigs.megahit.renamed.fa $cluster_dir/contigs.megahit.fa
    fi

    if [ ! -f $cluster_dir/contigs.megahit.fa.amb ]; then
        bwa index $cluster_dir/contigs.megahit.fa
    fi
    # if $reads2 is empty, then it is interleaved pair-end reads
    if [ -z $reads2 ]; then
        bwa mem -p -t $thread $cluster_dir/contigs.megahit.fa $reads1 | samtools sort -@ $thread -o $assembly_dir/contigs.megahit.bam     
    else
        bwa mem -t $thread $cluster_dir/contigs.megahit.fa $reads1 $reads2 | samtools sort -@ $thread -o $assembly_dir/contigs.megahit.bam
    fi

    jgi_summarize_bam_contig_depths --outputDepth $assembly_dir/contigs.megahit.depth $assembly_dir/contigs.megahit.bam &
    samtools sort -n -@ $thread $assembly_dir/contigs.megahit.bam -o $assembly_dir/contigs.megahit.name_sorted.bam &
    wait
    rm $assembly_dir/contigs.megahit.bam
fi