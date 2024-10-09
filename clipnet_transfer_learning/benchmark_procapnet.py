import glob
import os

import h5py
import numpy as np
import procapnet
import torch
from tangermeme.io import extract_loci

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

models = [procapnet.ProCapNet() for i in range(7)]
loaded_models = [
    models[i].load_state_dict(
        torch.load(glob.glob(f"../models/procapnet_k562/fold_{i}/*.torch")[0])
    )
    for i in range(7)
]

FOLDS = [
    ["chr1", "chr4"],
    ["chr2", "chr13", "chr16"],
    ["chr5", "chr6", "chr20", "chr21"],
    ["chr7", "chr8", "chr9"],
    ["chr10", "chr11", "chr12"],
    ["chr3", "chr14", "chr15", "chr17"],
    ["chr18", "chr19", "chr22", "chrX", "chrY"],
]

bed_path = "../../data/k562_procap_pairedPeak_autosomes.bed.gz"
genome_path = "../../data/hg38.fa"
pl_bw = "../../data/k562/k562_procap_pl.bigWig"
mn_bw = "../../data/k562/k562_procap_mn.bigWig"

loci_folds = [
    extract_loci(
        loci=bed_path,
        sequences=genome_path,
        signals=[pl_bw, mn_bw],
        chrom=fold,
    )
    for fold in FOLDS
]


# alt_pred = []
# for model in tqdm.tqdm(models):
#    alt_pred.append(predict(model, torch.tensor(alt_ohe).to(torch.float)))

alt_quantity = np.array([p[1][:, 0].numpy() for p in alt_pred])

with h5py.File("../data/mpra/k562_mpra_snps_2114_alt_procapnet_folds.h5", "w") as f:
    f.create_dataset("alt", data=np.array(alt_quantity), compression="gzip")
