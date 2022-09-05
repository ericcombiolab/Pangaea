import logging
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from rph_kmeans import RPHKMeans

from utils import run_cmd


def clustering_rph_kmeans(embedding, k):
    clt = RPHKMeans(n_init=20, n_clusters=k, verbose=0)
    clusters = clt.fit_predict(embedding)
    return clusters


def save_clustering_result(outdir, reads1, reads2, contigs, local_asm, athena, cutoff, clusters, barcodes, threads, script_path):
    snakemake = os.path.join(script_path, "bin", "Lathe", "Snakefile")
    extract_reads = os.path.join(script_path, "bin", "extract_reads")
    bin_assembly = os.path.join(script_path, "bin", "bin_assembly.sh")
    low_abd_reads = os.path.join(script_path, "bin", "low_abd_reads.sh")
    low_assembly = os.path.join(script_path, "bin", "low_assembly.sh")
    merge_asm = os.path.join(script_path, "bin", "merge_asm.sh")

    cluster_dir = os.path.join(outdir, "3.clustering")
    if not os.path.isdir(cluster_dir):
        run_cmd(["mkdir", cluster_dir])
        bc_npz = os.path.join(cluster_dir, "barcodes.npz")
        output_npz = os.path.join(cluster_dir, "clusters.npz")
        np.savez(bc_npz, barcodes)
        np.savez(output_npz, clusters)
        cluster2barcodes = defaultdict(list)
        for i in range(len(barcodes)):
            cluster2barcodes[clusters[i]].append(barcodes[i])
        output_path = os.path.join(cluster_dir, "clusters.tsv")
        with open(output_path, "w") as f:
            for c in cluster2barcodes:
                f.write("{}\t{}\n".format(c, ",".join(cluster2barcodes[c])))
        run_cmd([extract_reads, "-1", reads1, "-2", reads2, "-c", output_path, "-o", os.path.join(cluster_dir, "cluster")])
    
    assembly_dir = os.path.join(outdir, "4.assembly")
    cutoff = cutoff.split(',')
    if threads >= 150:
        threads = 150
    logging.info("mapping reads to contigs")
    run_cmd([bin_assembly, cluster_dir, assembly_dir, reads1, reads2, str(threads)])
    logging.info("obtaining reads mapped to low-abundance contigs")
    for cf in cutoff:
        run_cmd([low_abd_reads, cluster_dir, assembly_dir, str(cf)])
    logging.info("reassemble low-abundance contigs")
    executor = ThreadPoolExecutor(max_workers=5)
    for cf in cutoff:
        executor.submit(run_cmd, [low_assembly, assembly_dir, contigs, str(cf), str(threads)])
    executor.shutdown()
    logging.info("merge contigs with local assemblies")
    run_cmd([merge_asm, cluster_dir, assembly_dir, local_asm, athena, contigs, snakemake])
