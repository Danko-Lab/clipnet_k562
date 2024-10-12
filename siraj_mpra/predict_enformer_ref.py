"""
This script requires the enformer-pytorch package to be installed. (+ tangermeme)
"""

import h5py
import numpy as np
import pyfastx
import torch
import tqdm
from enformer_pytorch import from_pretrained
from tangermeme.predict import predict
from tangermeme.utils import one_hot_encode

# Load model

dnase_idx = [33, 34, 35, 121, 122, 123, 625]
cage_idx = [4828, 5111]


class K562Wrapper(torch.nn.Module):
    def __init__(self, model):
        super(K562Wrapper, self).__init__()
        self.model = model

    def forward(self, X):
        return self.model(X)["human"][:, :, dnase_idx + cage_idx]


enformer = K562Wrapper(from_pretrained("EleutherAI/enformer-official-rough"))

# Load sequences

fa = "../data/mpra/media-3_oligos_snps_cleaned_enformer_windows_holdouts.fa.gz"
seqs = [rec.seq.upper() for rec in tqdm.tqdm(pyfastx.Fasta(fa))]
ref = [one_hot_encode(seq) for seq in tqdm.tqdm(seqs)]
ref_tensor = torch.stack(ref).to(torch.float).swapaxes(1, 2)

# Predict

ref_prediction = predict(
    enformer, ref_tensor, device="cuda:0", batch_size=1, verbose=True
)
torch.cuda.empty_cache()

with h5py.File("../data/mpra/k562_mpra_snps_ref_enformer.h5", "w") as f:
    f.create_dataset("ref", data=np.array(ref_prediction), compression="gzip")
