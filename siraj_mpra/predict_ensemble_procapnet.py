import glob
import os
import sys

import h5py
import numpy as np
import procapnet
import torch
import tqdm

sys.path.append("../")
import utils

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

models = [
    procapnet.ProCapNet(
        glob.glob(f"../models/procapnet_k562/fold_{i}/*.torch")[0]
    ).cuda()
    for i in range(7)
]

ref_seqs = (
    utils.get_twohot_fasta_sequences(
        "../data/mpra/k562_mpra_snps_2114_ref.fa.gz"
    ).swapaxes(1, 2)
    / 2
)
ref_pred = []
for model in tqdm.tqdm(models):
    ref_pred.append(model.predict(torch.tensor(ref_seqs).to(torch.float).cuda()))

ref_quantity = np.mean([p[1] for p in ref_pred], axis=0)[:, 0]

alt_seqs = (
    utils.get_twohot_fasta_sequences(
        "../data/mpra/k562_mpra_snps_2114_alt.fa.gz"
    ).swapaxes(1, 2)
    / 2
)
alt_pred = []
for model in tqdm.tqdm(models):
    alt_pred.append(model.predict(torch.tensor(alt_seqs).to(torch.float).cuda()))

alt_quantity = np.mean([p[1] for p in alt_pred], axis=0)[:, 0]

with h5py.File("../data/mpra/k562_mpra_snps_2114_procapnet_ensemble.h5", "w") as f:
    f.create_dataset("ref", data=ref_quantity, compression="gzip")
    f.create_dataset("alt", data=alt_quantity, compression="gzip")
