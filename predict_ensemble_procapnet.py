import procapnet
import numpy as np
import h5py
import glob
import utils
import torch
import tqdm

paths = glob.glob("models/procapnet_k562/*.model")
models = [procapnet.Model(p).cuda() for p in paths]

ref_seqs = utils.get_twohot_fasta_sequences("data/mpra/k562_mpra_snps_ref.fa.gz").swapaxes(1, 2) / 2
alt_seqs = utils.get_twohot_fasta_sequences("data/mpra/k562_mpra_snps_alt.fa.gz").swapaxes(1, 2) / 2

ref_pred = []
for model in tqdm.tqdm(models):
    ref_pred.append(model.predict(torch.tensor(ref_seqs).to(torch.float).cuda()))

alt_pred = []
for model in tqdm.tqdm(models):
    alt_pred.append(model.predict(torch.tensor(alt_seqs).to(torch.float).cuda()))

ref_profile = np.mean([p[0] for p in ref_pred], axis=0)
ref_quantity = np.mean([p[1] for p in ref_pred], axis=0)[:, 0]

alt_profile = np.mean([p[0] for p in alt_pred], axis=0)
alt_quantity = np.mean([p[1] for p in alt_pred], axis=0)[:, 0]

with h5py.File("data/mpra/k562_mpra_snps_ref_procapnet.h5", "w") as f:
    f.create_dataset("profile", data=ref_profile, compression="gzip")
    f.create_dataset("quantity", data=ref_quantity, compression="gzip")

with h5py.File("data/mpra/k562_mpra_snps_alt_procapnet.h5", "w") as f:
    f.create_dataset("profile", data=alt_profile, compression="gzip")
    f.create_dataset("quantity", data=alt_quantity, compression="gzip")