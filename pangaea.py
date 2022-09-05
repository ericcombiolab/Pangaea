#!/usr/bin/env python
from ast import arg
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

from clustering import save_clustering_result
from data import Data
from feature import Feature
from models.VAENET import VAENET
from utils import CustomWeightedRandomSampler, init_all

def run(args, script_path):
    # init
    init_all(seed=2021, threads=args.threads, logfile="log", level=logging.INFO, outdir=args.output)
    logging.info("command: " + ' '.join(sys.argv))
    logging.info(args)

    barcodes, abundance, tnf = Feature(args.reads1, args.reads2, args.output, args.kmer, args.vector_size, args.window_size, args.min_length, args.threads, script_path).extract_features()

    # prepare data
    dataset = Data(barcodes, abundance, tnf)
    test_size = min(int(len(dataset)*0.7), 1000000)
    dataloader_train = DataLoader(dataset, batch_size=args.batch_size, num_workers=args.threads, shuffle=False, sampler=CustomWeightedRandomSampler(dataset.weights, num_samples=len(dataset)))
    dataloader_test = DataLoader(dataset, batch_size=args.batch_size, num_workers=args.threads, shuffle=False, sampler=CustomWeightedRandomSampler(dataset.weights, num_samples=test_size, replacement=False))
    dataloader_original = DataLoader(dataset, batch_size=args.batch_size, num_workers=args.threads, shuffle=True)

    # train
    model = VAENET(
        abd_dim=abundance.shape[1], tnf_dim=tnf.shape[1], latent_size=args.latent_dim, num_classes=args.clusters,
        epochs=args.epochs, cuda=args.use_cuda, num_gpus=args.num_gpus, lr=args.lr, dropout=args.dropout,
        alpha=args.weight_alpha, w_kl=args.weight_kl, weight_decay=args.weight_decay
    )
    labels, barcodes = model.train(dataloader_train, dataloader_test, dataloader_original, args.output, args.patience)
    save_clustering_result(args.output, args.reads1, args.reads2, args.spades, args.local_assembly, args.athena, args.low_abd_cut, labels, barcodes, args.threads, script_path)
    logging.info("program finished successfully")



def main():
    parser = argparse.ArgumentParser()
    # IO
    parser.add_argument("-1", "--reads1", required=True, default="", help="path to reads1 file (linked-reads)")
    parser.add_argument("-2", "--reads2", required=True, default="", help="path to reads2 file (linked-reads)")
    parser.add_argument("-o", "--output", required=True, help="output directory")
    # feature
    parser.add_argument("-l", "--min_length", type=int, default=2000, help="min barcode length (default 2000)")
    parser.add_argument("-k", "--kmer", type=int, default=15, help="kmer for abundance (default 15)")
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
    parser.add_argument("-sp", "--spades", type=str, required=True, help="path to original contigs")
    parser.add_argument("-lc", "--local_assembly", type=str, required=True, help="path to local assembly contigs")
    parser.add_argument("-at", "--athena", type=str, required=True, help="path to athena contigs")
    parser.add_argument("-lt", "--low_abd_cut", type=str, default="10,30,50,70,90", help="coverage for low abundance contigs")

    args = parser.parse_args()
    script_path = os.path.dirname(os.path.abspath(__file__))
    run(args, script_path)

    

if __name__ == "__main__":
    main()
