import procapnet
import numpy as np
import h5py
import glob
import utils
import torch

paths = glob.glob("models/procapnet_k562/*.model")
models = [procapnet.Model(p).cuda() for p in paths]

ref_seqs = utils.get_twohot_fasta_sequences("data/mpra/k562_mpra_snps_ref.fa.gz").swapaxes(1, 2) / 2
alt_seqs = utils.get_twohot_fasta_sequences("data/mpra/k562_mpra_snps_alt.fa.gz").swapaxes(1, 2) / 2


models[0].predict(torch.tensor(ref_seqs[:100]).to(torch.float).cuda())