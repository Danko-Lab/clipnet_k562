import glob
import os

import h5py
import torch
import tqdm
import numpy as np

import procapnet
import utils
from tangermeme.predict import predict

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

models = [
    procapnet.Model(glob.glob(f"models/procapnet_k562/fold_{i}/*.torch")[0]).cuda()
    for i in range(7)
]

ref_seqs = (
    utils.get_twohot_fasta_sequences(
        "data/mpra/k562_mpra_snps_2114_ref.fa.gz"
    ).swapaxes(1, 2)
    / 2
)
ref_pred = []
for model in tqdm.tqdm(models):
    ref_pred.append(predict(model, torch.tensor(ref_seqs).to(torch.float).cuda()))

ref_quantity = [p[1][:, 0] for p in ref_pred]


alt_seqs = (
    utils.get_twohot_fasta_sequences(
        "data/mpra/k562_mpra_snps_2114_alt.fa.gz"
    ).swapaxes(1, 2)
    / 2
)
alt_pred = []
for model in tqdm.tqdm(models):
    alt_pred.append(model.predict(torch.tensor(alt_seqs).to(torch.float).cuda()))

alt_quantity = [p[1][:, 0] for p in alt_pred]

with h5py.File("data/mpra/k562_mpra_snps_2114_procapnet_folds.h5", "w") as f:
    f.create_dataset("ref", data=np.array(ref_quantity), compression="gzip")
    f.create_dataset("alt", data=np.array(alt_quantity), compression="gzip")
