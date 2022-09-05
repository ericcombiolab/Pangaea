import logging
import numpy as np
import torch
from sklearn.preprocessing import normalize
from torch.utils.data import Dataset


class Data(Dataset):
    def __init__(self, barcodes, abd, tnf):
        super().__init__()
        self.bc = barcodes
        self.abd = abd
        self.tnf = tnf

        logging.info("calculate sampling weights")
        normalized_abd = normalize(self.abd, "l1")
        self.weights = np.array([(normalized_abd[i, ].max()) ** 2 for i in range(self.abd.shape[0])], dtype=np.float64)

        logging.info("normalize data")
        self.abd = normalized_abd.astype(np.float32)
        self.tnf = normalize(self.tnf, "l1").astype(np.float32)
        logging.info("preprocessing completed")


    def __len__(self):
        return self.abd.shape[0]

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        return {"abd": self.abd[idx, :], "tnf": self.tnf[idx, :], "bc": self.bc[idx]}
