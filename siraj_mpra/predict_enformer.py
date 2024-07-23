import sys

import h5py
import numpy as np
import torch
import tqdm
from enformer_pytorch import from_pretrained
from tangermeme.predict import predict
from tangermeme.utils import random_one_hot

sys.path.append("../")
import utils

# os.environ["CUDA_VISIBLE_DEVICES"] = "0"


class HumanWrapper(torch.nn.Module):
    def __init__(self, model):
        super(HumanWrapper, self).__init__()
        self.model = model

    def forward(self, X):
        return self.model(X)["human"]


enformer = HumanWrapper(from_pretrained("EleutherAI/enformer-official-rough"))
ref = random_one_hot((300_000, 4, 196_608), random_state=0).swapaxes(1, 2)

ref_prediction = predict(enformer, ref, device="cuda:0", batch_size=1, verbose=True)
torch.cuda.empty_cache()

with h5py.File("../data/mpra/k562_mpra_snps_ref_enformer.h5", "w") as f:
    f.create_dataset("ref", data=np.array(ref), compression="gzip")


alt_ohe = (
    utils.get_twohot_fasta_sequences(
        "../data/mpra/k562_mpra_snps_2114_alt.fa.gz"
    ).swapaxes(1, 2)
    / 2
)

alt_pred = []
for model in tqdm.tqdm(models):
    alt_pred.append(predict(model, torch.tensor(alt_ohe).to(torch.float)))

alt_quantity = np.array([p[1][:, 0].numpy() for p in alt_pred])

with h5py.File("../data/mpra/k562_mpra_snps_2114_alt_procapnet_folds.h5", "w") as f:
    f.create_dataset("alt", data=np.array(alt_quantity), compression="gzip")
