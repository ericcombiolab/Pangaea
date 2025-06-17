import logging
import os
import subprocess
import sys

import numpy as np
import torch
from torch.utils.data.sampler import WeightedRandomSampler


class CustomWeightedRandomSampler(WeightedRandomSampler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __iter__(self):
        rand_tensor = np.random.choice(
            range(0, len(self.weights)),
            size=self.num_samples,
            p=self.weights.numpy() / torch.sum(self.weights).numpy(),
            replace=self.replacement
        )
        rand_tensor = torch.from_numpy(rand_tensor)
        return iter(rand_tensor.tolist())


class EarlyStopping:
    def __init__(self, patience=7, delta=0, path='checkpoint.pt'):
        self.patience = patience
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta
        self.path = path

    def __call__(self, val_loss, model):
        score = -val_loss
        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score + self.delta:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        torch.save(model.state_dict(), self.path)
        self.val_loss_min = val_loss


def run_cmd_with_pipe(command, pipe_file=None):
    if not pipe_file:
        log_pipe = subprocess.DEVNULL
    else:
        # get pipe_file pipe
        log_pipe = open(pipe_file, "a")
    logging.info("command started: " + " ".join(command))
    # run command and check return code and save log to log_file
    ret = subprocess.run(command, stdout=log_pipe, stderr=log_pipe)
    if ret.returncode:
        logging.error("command failed: " + " ".join(command))
        sys.exit(1)
    logging.info("command completed: " + " ".join(command))

def run_cmd(command, log_file=None):
    if not log_file:
        log_pipe = subprocess.DEVNULL
    else:
        # get log_file pipe
        log_pipe = open(log_file, "a")
    logging.info("command started: " + " ".join(command))
    # run command and check return code and save log to log_file
    ret = subprocess.run(command, stdout= subprocess.DEVNULL, stderr=log_pipe)
    if ret.returncode:
        logging.error("command failed: " + " ".join(command))
        sys.exit(1)
    logging.info("command completed: " + " ".join(command))


def init_all(seed, threads, logfile, level, outdir):
    # random seed
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.set_num_threads(threads)
    # output directory
    if not os.path.isdir(outdir):
        subprocess.run(["mkdir", outdir])
    # root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    formatter = logging.Formatter("%(asctime)s (%(levelname)s): %(message)s", "%Y-%m-%d %H:%M:%S")
    handler1 = logging.FileHandler(os.path.join(outdir, logfile))
    handler2 = logging.StreamHandler()
    handler1.setLevel(level)
    handler2.setLevel(level)
    handler1.setFormatter(formatter)
    handler2.setFormatter(formatter)
    root_logger.addHandler(handler1)
    root_logger.addHandler(handler2)
    root_logger.info("program start up")
