import logging
import os

import numpy as np
import torch
import torch.nn as nn
from clustering import clustering_rph_kmeans
from torch.nn import functional as F
from utils import EarlyStopping, run_cmd


class VAENET:
    def __init__(self, abd_dim, tnf_dim, latent_size, num_classes, epochs, cuda, num_gpus, lr, dropout, alpha, w_kl, weight_decay):
        self.num_epochs = epochs
        self.num_classes = num_classes
        self.latent_size = latent_size
        self.input_size = abd_dim + tnf_dim
        self.learning_rate = lr
        self.weight_decay = weight_decay
        self.w_kl = w_kl*100/latent_size
        self.wa = alpha*100/np.log(abd_dim)
        self.wt = (1-alpha)*100/np.log(tnf_dim)
        self.cuda = cuda
        self.eps = 1e-9
        self.network = VaritionalAutoEncoder(abd_dim, tnf_dim, dropout=dropout)
        if self.cuda:
            self.network = self.network.cuda()
            self.network = nn.DataParallel(self.network, device_ids=range(num_gpus))
            self.network = self.network.module

    def train(self, train_loader, val_loader, dataloader, model_path, patience):
        if not os.path.isdir(model_path):
            # raise Exception and return
            raise Exception("model path not exist")
        train_model = os.path.join(model_path, "train_model.pk")
        early = EarlyStopping(patience=patience, delta=1e-6, path=train_model)
        if not os.path.exists(train_model):
            logging.info("train start")
            optimizer = torch.optim.Adam(self.network.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
            total_loss = []
            abd_loss = []
            tnf_loss = []
            kl_loss = []
            for epoch in range(1, self.num_epochs + 1):
                for batch, data in enumerate(train_loader):
                    self.network.train()
                    optimizer.zero_grad()
                    abd = data["abd"]
                    tnf = data["tnf"]
                    if self.cuda:
                        abd = abd.cuda()
                        tnf = tnf.cuda()
                    out_net = self.network(abd, tnf)
                    unlab_loss_dic = self.unlabeled_loss(out_net)
                    total = unlab_loss_dic['total']

                    total_loss.append(total.item())
                    abd_loss.append(unlab_loss_dic['abd_rec'].item())
                    tnf_loss.append(unlab_loss_dic['tnf_rec'].item())
                    kl_loss.append(unlab_loss_dic['kl_loss'].item())
                    total.backward()
                    optimizer.step()

                    if (batch + 1) % 100 == 0:
                        self.network.eval()
                        total_loss_val = []
                        barcodes = []
                        embedding = []
                        with torch.no_grad():
                            for data in val_loader:
                                abd = data["abd"]
                                tnf = data["tnf"]
                                if self.cuda:
                                    abd = abd.cuda()
                                    tnf = tnf.cuda()
                                out_net = self.network(abd, tnf)
                                unlab_loss_dic = self.unlabeled_loss(out_net)
                                total_loss_val.append(unlab_loss_dic['total'].item())
                                embedding.append(out_net["mu"].detach().cpu().numpy())
                                barcodes.extend(data["bc"])

                        logging.info(f"epoch {epoch}/{self.num_epochs} batch {batch+1}/{len(train_loader)}: train {np.average(total_loss):.8f} abd {np.average(abd_loss):.8f} tnf {np.average(tnf_loss):.8f} kl {np.average(kl_loss):.8f} | test {np.average(total_loss_val):.8f}")
                        early(np.average(total_loss_val), self.network)
                        total_loss = []
                        abd_loss = []
                        tnf_loss = []
                        kl_loss = []

                    if early.early_stop:
                        logging.info("early stop triggered")
                        break

                if early.early_stop:
                    logging.info("early stop triggered")
                    break
                self.network.eval()
                total_loss_val = []
                barcodes = []
                with torch.no_grad():
                    for data in val_loader:
                        abd = data["abd"]
                        tnf = data["tnf"]
                        if self.cuda:
                            abd = abd.cuda()
                            tnf = tnf.cuda()
                        out_net = self.network(abd, tnf)
                        unlab_loss_dic = self.unlabeled_loss(out_net)
                        total_loss_val.append(unlab_loss_dic['total'].item())
                        barcodes.extend(data["bc"])
                logging.info(f"epoch {epoch}/{self.num_epochs} batch {batch+1}/{len(train_loader)}: train {np.average(total_loss):.8f} abd {np.average(abd_loss):.8f} tnf {np.average(tnf_loss):.8f} kl {np.average(kl_loss):.8f} | test {np.average(total_loss_val):.8f}")
                total_loss = []
                abd_loss = []
                tnf_loss = []
                kl_loss = []

            if not os.path.exists(train_model):
                torch.save(self.network.state_dict(), train_model)
        else:
            logging.info("trainning model already saved")

        barcodes = []
        embedding = []
        latent_path = os.path.join(model_path, "latent.npz")
        barcodes_path = os.path.join(model_path, "barcodes.npz")
        if not os.path.exists(latent_path) or not os.path.exists(barcodes_path):
            self.network.load_state_dict(torch.load(train_model))
            self.network.eval()
            with torch.no_grad():
                for data in dataloader:
                    abd = data["abd"]
                    tnf = data["tnf"]
                    if self.cuda:
                        abd = abd.cuda()
                        tnf = tnf.cuda()
                    embedding.append(self.network.emebdding(abd, tnf).detach().cpu().numpy())
                    barcodes.extend(data["bc"])
            embedding = np.concatenate(embedding, axis=0)
            np.savez(barcodes_path, barcodes)
            np.savez(latent_path, embedding)
        else:
            logging.info("latent and barcodes already saved")
        # new a file show model finished
        with open(os.path.join(model_path, "model_finished"), "w") as f:
            f.write("model finished")
            
            #list(data.keys())
        #     embedding = np.load(latent_path)['arr_0']
        #     barcodes = np.load(barcodes_path)['arr_0']
        # logging.info("embedding.shape:")
        # logging.info(embedding.shape)

        # labels = clustering_rph_kmeans(embedding, self.num_classes)

        # return embedding, barcodes

    def unlabeled_loss(self, out_net):
        abd = out_net["abd"]
        tnf = out_net["tnf"]
        abd_rec = out_net["abd_rec"]
        tnf_rec = out_net["tnf_rec"]
        mu = out_net["mu"]
        logsigma = out_net["logsigma"]

        loss_abd = self.reconstruction_loss(abd, abd_rec)
        loss_tnf = self.reconstruction_loss(tnf, tnf_rec)
        loss_kl = -0.5 * (1 + logsigma - mu.pow(2) - logsigma.exp()).sum(dim=1).mean()
        loss_total = self.wa * loss_abd + self.wt * loss_tnf + self.w_kl * loss_kl

        loss_dic = {
            'total': loss_total,
            'abd_rec': loss_abd,
            'tnf_rec': loss_tnf,
            'kl_loss': loss_kl
        }
        return loss_dic

    def reconstruction_loss(self, real, predicted):
        loss = -(torch.log(predicted + self.eps)*real).sum(-1).mean()
        return loss


class VaritionalAutoEncoder(nn.Module):
    def __init__(
        self,
        input_abd_size,
        input_tnf_size,
        hidden_sizes=[512, 512],
        latent_size=32,
        dropout=0.2
    ):
        super().__init__()
        self.abd_size = input_abd_size
        self.tnf_size = input_tnf_size
        self.input_size = self.abd_size + self.tnf_size

        self.encoder = []
        for feature_size_in, feature_size_out in zip([self.input_size] + hidden_sizes, hidden_sizes):
            self.encoder.append(nn.Linear(feature_size_in, feature_size_out))
            self.encoder.append(nn.BatchNorm1d(feature_size_out))
            self.encoder.append(nn.LeakyReLU(True))
            self.encoder.append(nn.Dropout(dropout))
        self.encoder = nn.Sequential(*self.encoder)

        self.l_mu = nn.Linear(hidden_sizes[-1], latent_size)
        self.l_sigma = nn.Linear(hidden_sizes[-1], latent_size)
        self.softplus = nn.Softplus()

        self.decoder = []
        for feature_size_in, feature_size_out in zip([latent_size] + hidden_sizes[::-1], hidden_sizes[::-1]):
            self.decoder.append(nn.Linear(feature_size_in, feature_size_out))
            self.decoder.append(nn.BatchNorm1d(feature_size_out))
            self.decoder.append(nn.LeakyReLU(True))
            self.decoder.append(nn.Dropout(dropout))
        self.decoder = nn.Sequential(*self.decoder)
        self.output = nn.Linear(hidden_sizes[0], self.input_size)

    def calcu_latent(self, abd, tnf):
        input_tensor = torch.cat((abd, tnf), 1)
        hidden = self.encoder(input_tensor)
        mu = self.l_mu(hidden)
        logsigma = self.softplus(self.l_sigma(hidden))
        epsilon = torch.randn(mu.size(0), mu.size(1))
        epsilon.requires_grad = True
        latent = mu + epsilon * torch.exp(logsigma / 2)
        return mu, logsigma, latent

    def emebdding(self, abd, tnf):
        input_tensor = torch.cat((abd, tnf), 1)
        hidden = self.encoder(input_tensor)
        mu = self.l_mu(hidden)
        return mu

    def forward(self, abd, tnf):
        mu, logsigma, latent_tensor = self.calcu_latent(abd, tnf)
        output_tensor = self.output(self.decoder(latent_tensor))
        abd_rec = output_tensor.narrow(1, 0, self.abd_size)
        tnf_rec = output_tensor.narrow(1, self.abd_size, self.tnf_size)
        abd_rec = F.softmax(abd_rec, dim=1)
        tnf_rec = F.softmax(tnf_rec, dim=1)

        out_net = {}
        out_net["abd"] = abd
        out_net["tnf"] = tnf
        out_net["abd_rec"] = abd_rec
        out_net["tnf_rec"] = tnf_rec
        out_net["mu"] = mu
        out_net["logsigma"] = logsigma
        return out_net
