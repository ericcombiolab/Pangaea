import logging
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from sklearn.metrics import silhouette_score

import numpy as np
from rph_kmeans import RPHKMeans
from sklearn.cluster import KMeans

from utils import run_cmd, run_cmd_with_pipe
import warnings 

def clustering_rph_kmeans(embedding, k):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clt = RPHKMeans(n_init=20, n_clusters=k, verbose=0)
        clusters = clt.fit_predict(embedding)
        return clusters


def save_clustering_result(args, outdir, reads, clusters, barcodes, script_path):
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
        if isinstance(reads, str):
            run_cmd([extract_reads, "-i", reads, "-c", output_path, "-o", os.path.join(cluster_dir, "cluster")])
        elif isinstance(reads, list):
            run_cmd([extract_reads, "-1", reads[0], "-2", reads[1], "-c", output_path, "-o", os.path.join(cluster_dir, "cluster")])
    
    assembly_dir = os.path.join(outdir, "4.assembly")
    cutoff = args.low_abd_cut
    cutoff = cutoff.split(',')
    if args.threads >= 150:
        args.threads = 150
    logging.info("mapping reads to contigs")
    if isinstance(reads, str):
        run_cmd([bin_assembly, cluster_dir, assembly_dir, reads, str(args.threads)])
    elif isinstance(reads, list):
        run_cmd([bin_assembly, cluster_dir, assembly_dir, reads[0], reads[1], str(args.threads)])
    logging.info("obtaining reads mapped to low-abundance contigs")
    for cf in cutoff:
        run_cmd([low_abd_reads, cluster_dir, assembly_dir, str(cf)])
    logging.info("reassemble low-abundance contigs")
    executor = ThreadPoolExecutor(max_workers=5)
    for cf in cutoff:
        executor.submit(run_cmd, [low_assembly, assembly_dir, args.spades, str(cf), str(args.threads)])
    executor.shutdown()
    logging.info("merge contigs with local assemblies")
    run_cmd([merge_asm, cluster_dir, assembly_dir, args.local_assembly, args.athena, args.spades, snakemake])

def cluster_barcode_reads(args, model_path, cluster_path, script_path, range_step=(4, 100, 3), patience=7, delta=0.0001):    
    extract_reads = os.path.join(script_path, "bin", "extract_reads")
    calculate_diversity = os.path.join(script_path, "bin", "calculate_diversity.sh")
    output_npz = os.path.join(cluster_path, "clusters.npz")
    output_tsv = os.path.join(cluster_path, "clusters.tsv")
    embedding_path = os.path.join(model_path, "latent.npz")
    barcodes_path = os.path.join(model_path, "barcodes.npz")
    # check if latent.npz and barcodes.npz exists
    if not os.path.isfile(embedding_path) or not os.path.isfile(barcodes_path) or not os.path.isfile(extract_reads):
        raise FileNotFoundError(f"{embedding_path} or {barcodes_path} or {extract_reads} not found")
    if not os.path.isfile(output_tsv):
        try:
            # load embeding
            embedding = np.load(embedding_path)['arr_0']
            # load barcodes
            barcodes = np.load(barcodes_path)['arr_0']
            # clustering
            best_score = -1
            cluster2barcodes = defaultdict(list)

            if args.clusters:
                num_classes = args.clusters
            else:
                if args.reads1 and args.reads2:
                    run_cmd([calculate_diversity, args.reads1])
                elif args.interleaved_reads:
                    run_cmd([calculate_diversity, args.interleaved_reads])
                else:
                    logging.error("no reads provided")
                    raise FileNotFoundError("no reads provided")
                shannon_diversity = np.loadtxt("metaphlan_tmp/diversity_analysis/profiles_table_shannon.txt")
                num_classes = int(8 * shannon_diversity)
                logging.info(f"estimated num_classes: {num_classes}")
                run_cmd(["rm", "-rf", "metaphlan_tmp"])
            clusters = clustering_rph_kmeans(embedding, num_classes)
            np.savez(output_npz, clusters)
            for i in range(len(barcodes)):
                cluster2barcodes[clusters[i]].append(barcodes[i])
            logging.info("saving clustering tsv")
            
            with open(output_tsv, "w") as tsv:
                for cluster_id in cluster2barcodes:
                    tsv.write("{}\t{}\n".format(cluster_id, ",".join(cluster2barcodes[cluster_id])))
        except:
            logging.error("clustering failed")
            raise
    else:
        logging.info("existing clustering result found")
    
    if args.reads1 and args.reads2:
        run_cmd([extract_reads, "-1", args.reads1, "-2", args.reads2, "-c", output_tsv, "-o", os.path.join(cluster_path, "cluster")])
    elif args.interleaved_reads:
        run_cmd([extract_reads, "-i", args.interleaved_reads, "-c", output_tsv, "-o", os.path.join(cluster_path, "cluster")])
    
    else:
        logging.error("no reads provided")
        raise FileNotFoundError("no reads provided")
    with open(os.path.join(cluster_path, "clustering_finished"), "w") as f:
        # new a file show finished
        f.write("finished")


def final_assemble(args, cluster_path, assembly_path, script_path):
    snakemake = os.path.join(script_path, "bin", "Lathe", "Snakefile")
    bin_assembly = os.path.join(script_path, "bin", "bin_assembly.sh")
    low_abd_reads = os.path.join(script_path, "bin", "low_abd_reads.sh")
    low_assembly = os.path.join(script_path, "bin", "low_assembly.sh")
    merge_asm = os.path.join(script_path, "bin", "merge_asm.sh")
    threads = args.threads
    cutoff = args.low_abd_cut
    cutoff = cutoff.split(',')

    if threads >= 150:
        threads = 150
    logging.info("mapping reads to contigs")
    if args.reads1 and args.reads2:
        run_cmd([bin_assembly, cluster_path, assembly_path, str(threads), args.reads1, args.reads2])
    elif args.interleaved_reads:
        run_cmd([bin_assembly, cluster_path, assembly_path, str(threads), args.interleaved_reads])
    else :
        logging.error("no linked reads or interleaved reads provided")
        raise FileNotFoundError("no reads provided")
    logging.info("obtaining reads mapped to low-abundance contigs")
    for cf in cutoff:
        run_cmd([low_abd_reads, cluster_path, assembly_path, str(cf)])
    logging.info("reassemble low-abundance contigs")
    executor = ThreadPoolExecutor(max_workers=5)
    for cf in cutoff:
        executor.submit(run_cmd, [low_assembly, assembly_path, args.spades, str(cf), str(threads), args.low_assembler])
    executor.shutdown()
    logging.info("merge args.spades/(seed) contigs with local assemblies")
    run_cmd([merge_asm, cluster_path, assembly_path, args.local_assembly, args.athena, args.spades, args.low_assembler, snakemake])
    # new a file show finished
    with open(os.path.join(assembly_path, "assembly_finished"), "w") as f:
        f.write("finished")