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

fa = "/workdir/ayh8/data/k562_siraj_mpra/media-3_oligos_snps_cleaned_enformer_windows_holdouts_all.fa.gz"
seqs = {rec.name: rec.seq.upper() for rec in tqdm.tqdm(pyfastx.Fasta(fa))}
padded_seqs = {}
for name, seq in tqdm.tqdm(
    seqs.items(), total=len(seqs), desc="Padding sequences with Ns"
):
    if len(seq) == 196_608:
        padded_seqs[name] = seq
    else:
        padding = 196_608 - len(seq)
        start = name.split(":")[-1].split("-")[0]
        if start == "0":
            print(f"Padding start of {name} with {padding} Ns")
            seq = "N" * padding + seq
        else:
            print(f"Padding end of {name} with {padding} Ns")
            seq = seq + "N" * padding
        padded_seqs[name] = seq

alt_seqs = [
    seq[: len(seq) // 2 - 1] + name.split(":")[-1] + seq.upper()[len(seq) // 2 :]
    for name, seq in tqdm.tqdm(padded_seqs.items())
]
alt = [one_hot_encode(seq) for seq in tqdm.tqdm(alt_seqs, total=len(alt_seqs))]
alt_tensor = torch.stack(alt).to(torch.float).swapaxes(1, 2)

# Predict

alt_prediction = predict(
    enformer, alt_tensor, device="cuda:1", batch_size=1, verbose=True
)

torch.cuda.empty_cache()

with h5py.File(
    "/workdir/ayh8/data/k562_siraj_mpra/k562_mpra_snps_alt_enformer_all.h5", "w"
) as f:
    f.create_dataset("alt", data=np.array(alt_prediction), compression="gzip")
