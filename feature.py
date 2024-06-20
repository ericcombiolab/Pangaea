import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from utils import run_cmd, run_cmd_with_pipe


class Feature:
    def __init__(self, args,  script_path):
        # self.reads = args.reads
        self.args = args
        self.tnf_k = str(args.tnf_kmer)
        self.ws = args.window_size
        self.vs = args.vector_size
        self.kmer = args.kmer
        self.minl = args.min_length
        self.threads = args.threads
        self.count_kmer = os.path.join(script_path, "bin", "count_kmer")
        self.count_tnf = os.path.join(script_path, "bin", "count_tnf")
        self.jellyfish = "jellyfish"
        self.feature_dir = os.path.join(args.output, "1.features")
        if not os.path.isdir(self.feature_dir):
            run_cmd(["mkdir", self.feature_dir])
    
    def extract_features(self):
        executor = ThreadPoolExecutor(max_workers=3)
        abundance = executor.submit(self.run_jellyfish)
        tnf = executor.submit(self.calcu_tnf)
        readnames1, abundance = abundance.result()
        readnames2, tnf = tnf.result()
        executor.shutdown()
        assert (readnames1 == readnames2).all()
        # new a file show feature finished
        with open(os.path.join(self.feature_dir, "feature_finished"), "w") as f:
            f.write("feature finished")
        return readnames1, abundance, tnf
    
    def load_features(self):
        out_pkl = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.v{self.vs}.w{self.ws}.m{self.minl}.pkl")
        out_freq = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.v{self.vs}.w{self.ws}.m{self.minl}.gz")
        out_tnf = os.path.join(self.feature_dir, f"tnf.m{self.minl}.pkl")
        try:
            tnf = pd.read_pickle(out_tnf)
            readnames = tnf[0].to_numpy()
            tnf = tnf.drop(columns=0).to_numpy()
            logging.info(f"tnf shape {tnf.shape}")
        except:
            raise Exception(out_tnf, " file not found")

        if os.path.isfile(out_pkl):
            logging.info("load features from pickle file "+out_pkl)
            df = pd.read_pickle(out_pkl)
            readnames = df[0].to_numpy()
            abundance = df.drop(columns=0).to_numpy()
        elif os.path.isfile(out_freq):
            logging.info(out_pkl+" not found, load features from "+out_freq)
            readnames = pd.read_csv(out_freq, sep="\t", header=None, usecols=[0]).values.ravel()
            abundance = pd.read_csv(out_freq, sep="\t", header=None, usecols=range(1, self.vs + 1)).values
        else:
            raise Exception(out_freq, " file not found")
        logging.info(f"abundance shape {abundance.shape}")
        return readnames, abundance, tnf

    def run_jellyfish(self):
        out_pkl = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.v{self.vs}.w{self.ws}.m{self.minl}.pkl")
        out_freq = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.v{self.vs}.w{self.ws}.m{self.minl}.gz")
        out_dump = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.dump")
        out_count = os.path.join(self.feature_dir, f"abundance.k{self.kmer}.count")
        if not os.path.isfile(out_count):
            logging.info("caculate abundance : "+out_count )
            if self.args.reads1 and self.args.reads2:
                command = ["pigz", "-dc", self.args.reads1, self.args.reads2, "|", "jellyfish", "count", "-t", str(self.threads), "-C", "-m", str(self.kmer), "-s", "5G", "-o", out_count, "--min-qual-char=?", "/dev/fd/0"]
                logging.info("command started: " + " ".join(command))
                pipe = subprocess.Popen(["pigz", "-dc", self.args.reads1, self.args.reads2], stdout=subprocess.PIPE)
                subprocess.check_output(["jellyfish", "count", "-t", str(self.threads), "-C", "-m", str(self.kmer), "-s", "5G", "-o", out_count, "--min-qual-char=?", "/dev/fd/0", ], stdin=pipe.stdout)
                pipe.communicate()
                logging.info("command completed: " + " ".join(command))
                run_cmd(["jellyfish", "dump", "-c", "-t", out_count, "-o", out_dump])
            elif self.args.interleaved_reads:
                if self.args.interleaved_reads.endswith(".gz"):
                    out_fq = os.path.join(self.feature_dir, "interleaved_reads.fq")
                    run_cmd_with_pipe(["pigz", "-dc", self.args.interleaved_reads], out_fq)
                else:
                    out_fq = self.args.interleaved_reads
                run_cmd([self.jellyfish, "count", out_fq, "-t", str(self.threads), "-C", "-m", str(self.kmer), "-s", "5G", "-o", out_count])
                if self.args.interleaved_reads.endswith(".gz"):
                    run_cmd(["rm", out_fq])
                
            else:
                raise ValueError("reads must be specified")
            
        if not os.path.isfile(out_dump):
            logging.info("dump abundance : "+out_dump )
            run_cmd([self.jellyfish, "dump", "-c", "-t", out_count, "-o", out_dump])
        if not os.path.isfile(out_freq):
            logging.info("caculate frequency : "+out_freq )
            if self.args.reads1 and self.args.reads2:
                run_cmd([self.count_kmer, "-1", self.args.reads1,"-2", self.args.reads2, "-t", str(self.threads), "-g", out_dump, "-k", str(self.kmer), "-l", str(self.minl), "-w", str(self.ws), "-v", str(self.vs), "-o", out_freq])
            elif self.args.interleaved_reads:
                run_cmd([self.count_kmer, "-i", str(self.args.interleaved_reads), "-t", str(self.threads), "-g", out_dump, "-k", str(self.kmer), "-l", str(self.minl), "-w", str(self.ws), "-v", str(self.vs), "-o", out_freq])
            else:
                raise ValueError("reads must be specified")
                    
        if not os.path.isfile(out_pkl):
            logging.info("save abundance")
            abundance = pd.read_csv(out_freq, header=None)
            abundance.to_pickle(out_pkl)
        else:
            logging.info("load abundance")
            abundance = pd.read_pickle(out_pkl)
        readnames = abundance[0].to_numpy()
        abundance = abundance.drop(columns=0).to_numpy()
        logging.info(f"abundance shape {abundance.shape}")
        return readnames, abundance

    def calcu_tnf(self):
        out_freq = os.path.join(self.feature_dir, f"tnf.m{self.minl}.gz")
        out_pkl = os.path.join(self.feature_dir, f"tnf.m{self.minl}.pkl")
        
        if not os.path.isfile(out_freq):
            if self.args.reads1 and self.args.reads2:
                run_cmd([self.count_tnf, "-1", self.args.reads1, "-2", self.args.reads2, "-k", self.tnf_k, "-t", str(self.threads), "-l", str(self.minl), "-o", out_freq])
            elif self.args.interleaved_reads:
                run_cmd([self.count_tnf, "-i", self.args.interleaved_reads, "-k", self.tnf_k, "-t", str(self.threads), "-l", str(self.minl), "-o", out_freq])
            else:
                raise ValueError("reads must be specified")
            
            logging.info("caculate tnf")
        if not os.path.isfile(out_pkl):
            tnf = pd.read_csv(out_freq, header=None)
            tnf.to_pickle(out_pkl)
        else:
            logging.info("load tnf")
            tnf = pd.read_pickle(out_pkl)
        readnames = tnf[0].to_numpy()
        tnf = tnf.drop(columns=0).to_numpy()
        logging.info(f"tnf shape {tnf.shape}")
        return readnames, tnf