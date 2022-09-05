import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from utils import run_cmd


class Feature:
    def __init__(self, reads1, reads2, output, kmer, vector_size, window_size, min_length, threads, script_path):
        self.r1 = reads1
        self.r2 = reads2
        self.ws = window_size
        self.vs = vector_size
        self.kmer = kmer
        self.minl = min_length
        self.threads = threads
        self.count_kmer = os.path.join(script_path, "bin", "count_kmer")
        self.count_tnf = os.path.join(script_path, "bin", "count_tnf")
        self.feature_dir = os.path.join(output, "1.features")
        if not os.path.isdir(self.feature_dir):
            run_cmd(["mkdir", self.feature_dir])
    
    def extract_features(self):
        executor = ThreadPoolExecutor(max_workers=2)
        abundance = executor.submit(self.run_jellyfish)
        tnf = executor.submit(self.calcu_tnf)
        barcodes1, abundance = abundance.result()
        barcodes2, tnf = tnf.result()
        executor.shutdown()
        assert (barcodes1 == barcodes2).all()
        return barcodes1, abundance, tnf

    def run_jellyfish(self):
        out_pkl = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.v{self.vs}.w{self.ws}.m{self.minl}.pkl")
        out_freq = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.v{self.vs}.w{self.ws}.m{self.minl}.gz")
        out_dump = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.dump")
        out_count = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.count")
        if not os.path.isfile(out_pkl):
            if not os.path.isfile(out_dump):
                command = ["pigz", "-dc", self.r1, self.r2, "|", "jellyfish", "count", "-t", str(self.threads), "-C", "-m", str(self.kmer), "-s", "5G", "-o", out_count, "--min-qual-char=?", "/dev/fd/0"]
                logging.info("command started: " + " ".join(command))
                pipe = subprocess.Popen(["pigz", "-dc", self.r1, self.r2], stdout=subprocess.PIPE)
                subprocess.check_output(["jellyfish", "count", "-t", str(self.threads), "-C", "-m", str(self.kmer), "-s", "5G", "-o", out_count, "--min-qual-char=?", "/dev/fd/0", ], stdin=pipe.stdout)
                pipe.communicate()
                logging.info("command completed: " + " ".join(command))
                run_cmd(["jellyfish", "dump", "-c", "-t", out_count, "-o", out_dump])
            if not os.path.isfile(out_freq):
                run_cmd([self.count_kmer, "-1", self.r1, "-2", self.r2, "-t", str(self.threads), "-g", out_dump, "-k", str(self.kmer), "-l", str(self.minl), "-w", str(self.ws), "-v", str(self.vs), "-o", out_freq])
            logging.info("load abundance")
            abundance = pd.read_csv(out_freq, header=None)
            abundance.to_pickle(out_pkl)
        else:
            logging.info("load abundance")
            abundance = pd.read_pickle(out_pkl)
        barcodes = abundance[0].to_numpy()
        abundance = abundance.drop(columns=0).to_numpy()
        logging.info(f"abundance shape {abundance.shape}")
        return barcodes, abundance

    def calcu_tnf(self):
        out_freq = os.path.join(self.feature_dir, f"tnf.m{self.minl}.gz")
        out_pkl = os.path.join(self.feature_dir, f"tnf.m{self.minl}.pkl")
        if not os.path.isfile(out_pkl):
            if not os.path.isfile(out_freq):
                run_cmd([self.count_tnf, "-1", self.r1, "-2", self.r2, "-k", str(4), "-t", str(self.threads), "-l", str(self.minl), "-o", out_freq])
            logging.info("load tnf")
            tnf = pd.read_csv(out_freq, header=None)
            tnf.to_pickle(out_pkl)
        else:
            logging.info("load tnf")
            tnf = pd.read_pickle(out_pkl)
        barcodes = tnf[0].to_numpy()
        tnf = tnf.drop(columns=0).to_numpy()
        logging.info(f"tnf shape {tnf.shape}")
        return barcodes, tnf
