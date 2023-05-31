#!/usr/bin/env python
import os
import shutil
import subprocess
import sys

import pysam
script_path = os.path.dirname(os.path.abspath(__file__))

def concat_files(input_list, output_path, appendnl=True):
    with open(output_path, "w") as outf:
        for path in input_list:
            with open(path) as inf:
                nonempty = False
                for line in inf:
                    nonempty = True
                    outf.write(line)
                if appendnl and nonempty and not line.endswith("\n"):
                    print("WARNING last char not newline", path)
                    outf.write("\n")


def get_fasta_sizes(fa_path):
    fasta = pysam.FastaFile(fa_path)
    ctg_size_map = {}
    for ctg in fasta.references:
        size = fasta.get_reference_length(ctg)
        ctg_size_map[ctg] = size
    return ctg_size_map

#
def filter_inputs(mergedbam, local, filtfa):
    ctg_size_map = get_fasta_sizes(local)
    fhandle = pysam.Samfile(mergedbam)
    full_ctgs = set()
    for read in fhandle:
        if read.is_unmapped:
            continue
        if read.query_alignment_length + 1000 >= ctg_size_map[read.qname]:
            full_ctgs.add(read.qname)
    fhandle.close()

    new_ctgs = set(ctg_size_map.keys()) - full_ctgs
    fasta = pysam.FastaFile(local)
    with open(filtfa, "w") as fout:
        for ctg in new_ctgs:
            seq = str(fasta.fetch(ctg).upper())
            fout.write(">{}\n".format(ctg))
            fout.write("{}\n".format(seq))
    num_orig = len(ctg_size_map)
    num_filt = len(new_ctgs)
    return num_orig, num_filt


def run(seeds, local, outdir):
    if not os.path.isdir(outdir):
        os.system("mkdir " + out)
    filtfa = os.path.join(outdir, "pre-flye-input-contigs.filt.fa")
    dup_5_seedsfa = os.path.join(outdir, "seed-contigs.fa")
    mergedfa = os.path.join(outdir, "flye-input-contigs.fa")

    # align local(athena_local+lowabd+binning) to templates(spades(aka:seeds_ctgs))
    aling_bam = os.path.join(outdir, "align-inputs.bam")
    if not os.path.isfile(aling_bam):
        with open(os.devnull, "w") as devnull:
            if not os.path.isfile(seeds + ".amb"):
                os.system("bwa index " + seeds)
            cmd = "bwa mem -t {} {} {} | samtools view -bS - | samtools sort -o {} - ".format(100, seeds, local, aling_bam)
            pp = subprocess.Popen(cmd, shell=True, stdout=devnull, stderr=devnull)
            retcode = pp.wait()
            assert (retcode == 0), "bwa alignment of unfiltered subassembly contigs {} to {} failed".format(local, seeds)
            cmd = "samtools index {}".format(aling_bam)
            pp = subprocess.Popen(cmd, shell=True, stdout=devnull, stderr=devnull)
            retcode = pp.wait()
            assert retcode == 0, "indexing of {} failed".format(aling_bam)

    
    if not os.path.isfile(mergedfa):
        # filter out contigs that are fully covered by the templates
        filter_inputs(aling_bam, local, filtfa)
        os.system("seqtk seq -L 1000 {} > {}".format(seeds, dup_5_seedsfa))
        for _ in range(5):
            os.system("seqtk seq -L 1000 {} >> {}".format(seeds, dup_5_seedsfa))
        # merge the filtered contigs with the templates*5
        concat_files([filtfa, dup_5_seedsfa], mergedfa)
        os.system("{} {} contig_ > {}".format(os.path.join(script_path, "parse_header"), mergedfa, mergedfa+'.tmp'))
        os.system("mv {} {}".format(mergedfa+'.tmp', mergedfa))

    flye0_path = os.path.join(outdir, "flye-asm-1")
    flye_contigs_path = os.path.join(flye0_path, "assembly.fasta")
    cmd = "flye --meta --subassemblies {} --out-dir {} --threads {} --min-overlap 1000".format(mergedfa, flye0_path, 128)
    if not os.path.isfile(flye_contigs_path):
        subprocess.check_call(cmd, shell=True)

    final_fa_path = os.path.join(outdir, "final.asm.fa")
    shutil.copy(flye_contigs_path, final_fa_path)


templates = sys.argv[1]
local = sys.argv[2]
out = sys.argv[3]
# templates: spades(aka:seeds_ctgs)
# local: athena_local+lowabd+binning
run(templates, local, out)
