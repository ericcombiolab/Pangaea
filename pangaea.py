#!/usr/bin/env python
import os
DEFAULT_PROCESSES = os.cpu_count()
os.environ["MKL_NUM_THREADS"] = str(DEFAULT_PROCESSES)
os.environ["NUMEXPR_NUM_THREADS"] = str(DEFAULT_PROCESSES)
os.environ["OMP_NUM_THREADS"] = str(DEFAULT_PROCESSES)
os.environ["NUMEXPR_MAX_THREADS"] = str(DEFAULT_PROCESSES)

import sys
import argparse
import logging

from torch.utils.data import DataLoader
import numpy as np

from clustering import cluster_barcode_reads, final_assemble, final_assemble_long
from data import Data
from feature import Feature
from models.VAENET import VAENET
from utils import CustomWeightedRandomSampler, init_all, run_cmd

# check if step are finished
def check_steps_finish(args, step):
    if step == "1":
        return os.path.exists(os.path.join(args.output, "1.features")) and os.path.exists(os.path.join(args.output, "1.features", "feature_finished"))
    elif step == "2":
        if args.model == "vae":
            model_path = os.path.join(args.output, "2.vae")
        return  os.path.exists(model_path) and os.path.exists(os.path.join(model_path, "model_finished"))
    elif step == "3":
        return os.path.exists(os.path.join(args.output, "3.clustering")) and os.path.exists(os.path.join(args.output, "3.clustering", "clustering_finished"))
    elif step == "4":
        return os.path.exists(os.path.join(args.output, "4.assembly")) and os.path.exists(os.path.join(args.output, "4.assembly", "assemble_finished"))
    else:
        return False

# check if steps are required
def check_steps_required(args_steps, step):
    if step == "1":
        return "1" in args_steps
    elif step == "2":
        return "2" in args_steps
    elif step == "3":
        return "3" in args_steps
    elif step == "4":
        return "4" in args_steps
    else:
        return False

def run(args, script_path):
    # init
    init_all(seed=2021, threads=args.threads, logfile="log", level=logging.INFO, outdir=args.output)
    logging.info("command: " + ' '.join(sys.argv))
    logging.info(args)

    if  args.model == "vae":
        model_path = os.path.join(args.output, "2.vae")
    cluster_path = os.path.join(args.output, "3.clustering")
    assembly_path = os.path.join(args.output, "4.assembly")
    # step 1: feature extraction
    # new numpy array: read_specify, abundance, tnf
    read_specify, abundance, tnf = None, None, None
    if not check_steps_required(args.steps, "1"):
        logging.info("skip step 1: feature extraction")
    elif check_steps_finish(args, "1"):
        logging.info("step 1: feature extraction finished")
    else:
        # extract features
        if (args.reads1 and args.reads2) or args.interleaved_reads or args.long_reads:
            read_specify, abundance, tnf = Feature( args, script_path).extract_features()
        else:
            print("Please provide one or two input file(s):-1 and -2 for pair-end linked reads; -lr as long reads; -i for interleaved linked reads.")
            exit()

    # step 2: training
    if not check_steps_required(args.steps, "2"):
        logging.info("skip step 2: training")
    elif check_steps_finish(args, "2"):
        logging.info("step 2: training finished")
    else:
        # prepare data
        # if read_specify, abundance or tnf are not numpy array
        if not isinstance(read_specify, np.ndarray) or not isinstance(abundance, np.ndarray) or not isinstance(tnf, np.ndarray):
            read_specify, abundance, tnf = Feature(args, script_path).load_features()
        dataset = Data(read_specify, abundance, tnf) # read_specify could be barcodes(for linked reads) or readnames(for long reads)
        test_size = min(int(len(dataset)*0.7), 1000000)
        dataloader_train = DataLoader(dataset, batch_size=args.batch_size, num_workers=args.threads, shuffle=False, sampler=CustomWeightedRandomSampler(dataset.weights, num_samples=len(dataset)))
        dataloader_test = DataLoader(dataset, batch_size=args.batch_size, num_workers=args.threads, shuffle=False, sampler=CustomWeightedRandomSampler(dataset.weights, num_samples=test_size, replacement=False))
        dataloader_original = DataLoader(dataset, batch_size=args.batch_size, num_workers=args.threads, shuffle=True)
        if args.model == "vae":
            # check if vae exists
            if not os.path.exists(model_path):
                run_cmd(["mkdir", model_path])
            vae = VAENET(
                abd_dim=abundance.shape[1], tnf_dim=tnf.shape[1], latent_size=args.latent_dim, num_classes=args.clusters,
                epochs=args.epochs, cuda=args.use_cuda, num_gpus=args.num_gpus, lr=args.lr, dropout=args.dropout,
                alpha=args.weight_alpha, w_kl=args.weight_kl, weight_decay=args.weight_decay
            )
            vae.train(dataloader_train, dataloader_test, dataloader_original, model_path, args.patience)

    # step 3: clustering
    if not check_steps_required(args.steps, "3"):
        logging.info("skip step 3: clustering")
    elif check_steps_finish(args, "3"):
        logging.info("step 3: clustering finished")
    else:
        logging.info("start clustering")
        if not os.path.exists(cluster_path):
            run_cmd(["mkdir", cluster_path])
        if (args.reads1 and args.reads2) or args.long_reads or args.interleaved_reads:
            cluster_barcode_reads(args, model_path, cluster_path, script_path)
        else:
            logging.info("Please provide one or two input file(s):-1 and -2 for pair-end linked reads; -lr as long reads; -i for interleaved linked reads.")
            exit()
    
    # step 4: assembly
    if not check_steps_required(args.steps, "4"):
        logging.info("skip step 4: assembly")
    elif check_steps_finish(args, "4"):
        logging.info("step 4: assembly finished")
    else:
        logging.info("start assembly")
        if (args.reads1 and args.reads2) or args.interleaved_reads:
            final_assemble(args, cluster_path, assembly_path, script_path )
        elif args.long_reads:
            final_assemble_long(args, cluster_path, assembly_path, script_path )
    logging.info("program finished successfully")



def main():
    parser = argparse.ArgumentParser()
    # IO
    parser.add_argument("-1", "--reads1", default="", help="path to reads1 file (linked-reads)")
    parser.add_argument("-2", "--reads2", default="", help="path to reads2 file (linked-reads)")
    parser.add_argument("-lreads", "--long_reads", default="", help="path to reads file (long-reads)")
    #long reads type (pacbio_raw or nanopore_raw or pacbio_corrected or nanopore_corrected)
    parser.add_argument("-lrtype", "--long_reads_type", default="pacbio_raw", choices=["pacbio_raw", "nanopore_raw", "pacbio_corrected", "nanopore_corrected"], help="long reads type (default pacbio_raw)")
    parser.add_argument("-i", "--interleaved_reads", default="", help="path to reads file (long-reads)")

    parser.add_argument("-o", "--output", required=True, help="output directory")
    # feature
    parser.add_argument("-l", "--min_length", type=int, default=2000, help="min barcode length (default 2000)")
    parser.add_argument("-k", "--kmer", type=int, default=15, help="kmer for abundance (default 15)")
    parser.add_argument("-tnf_k", "--tnf_kmer", type=int, default=4, help="kmer for TNF (default 4, long reads should use 3)")
    parser.add_argument("-s", "--window_size", type=int, default=10, help="window size for abundance (default 10)")
    parser.add_argument("-v", "--vector_size", type=int, default=400, help="vector size for abundance (default 400)")
    # model
    parser.add_argument("-r", "--lr", type=float, default=0.005, help="learning rate (default 0.005)")
    parser.add_argument("-w", "--weight_decay", type=float, default=0.0001, help="weight decay (default 0.0001)")
    parser.add_argument("-e", "--epochs", type=int, default=100, help="number of epochs (default 100)")
    parser.add_argument("-b", "--batch_size", type=int, default=2048, help="batch size (defult 2048)")
    parser.add_argument("-d", "--dropout", type=float, default=0.2, help="dropout (default 0.2)")
    parser.add_argument("-p", "--patience", type=int, default=20, help="early stop patience (default 20)")
    parser.add_argument("-wa", "--weight_alpha", type=float, default=0.1, help="training weight for abundance and tnf (default 0.1)")
    parser.add_argument("-wk", "--weight_kl", type=float, default=0.015, help="training weight for KL (default 0.015)")
    parser.add_argument("-ld", "--latent_dim", type=int, default=32, help="latent dimension (default 32)")
    # others
    parser.add_argument("-c", "--clusters", type=int, required=True, help="number of clusters")
    parser.add_argument("-t", "--threads", type=int, default=100, help="number of threads (default 100)")
    parser.add_argument("-g", "--use_cuda", type=bool, default=False, help="use cuda (default False)")
    parser.add_argument("-n", "--num_gpus", type=int, default=1, help="use gpu in parallel (if use cuda)")
    parser.add_argument("-sp", "--spades", type=str, help="path to original contigs")
    parser.add_argument("-lc", "--local_assembly", type=str, help="path to local assembly contigs")
    parser.add_argument("-at", "--athena", type=str,  help="path to athena contigs")
    parser.add_argument("-lt", "--low_abd_cut", type=str, default="10,30,50,70,90", help="coverage for low abundance contigs")

    parser.add_argument("-md", "--model", type=str, default="vae", help="model ( vae)")
    parser.add_argument("-wx", "--weight_auxiliary", type=float, default=0.1, help="training weight for auxiliary (default 0.1)")
    parser.add_argument("-ls", "--loss_type", type=str, default="ce", help="reconstruction loss type (default ce)")
    parser.add_argument("-cf", "--confidence", type=float, default=0.85, help="clustering confidence")

    #step control
    parser.add_argument("-st", "--steps", type=str, default="1,2,3,4", help="steps to run (default 1:feature extraction, 2:vae trainning, 3:clutsering, 4:sub-assembly and final assembly)")

    args = parser.parse_args()
    script_path = os.path.dirname(os.path.abspath(__file__))
    run(args, script_path)

    

if __name__ == "__main__":
    main()
