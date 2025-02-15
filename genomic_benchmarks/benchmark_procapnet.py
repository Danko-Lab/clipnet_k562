"""
This script benchmarks the ProCapNet model on the K562 ProCap dataset.
It requires the ProCapNet model to be trained and saved in the models directory.
This also requires the ProCapNet environment to be installed.
"""

import glob
import os

import numpy as np
import pandas as pd
import procapnet
import torch
from scipy.stats import pearsonr
from tangermeme.io import extract_loci
from tangermeme.predict import predict

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

bed_path = "../../data/k562/k562_procap_pairedPeak_autosomes.bed.gz"
genome_path = "../../data/hg38.fa"
pl_bw = "../../data/k562/k562_procap_pl.bigWig"
mn_bw = "../../data/k562/k562_procap_mn.bigWig"

loci_folds = [
    extract_loci(
        loci=bed_path,
        sequences=genome_path,
        signals=[pl_bw, mn_bw],
        chroms=fold,
        verbose=True,
    )
    for fold in FOLDS
]

predictions = [
    predict(model, data[0].type(torch.float32), verbose=True)
    for model, data in zip(models, loci_folds)
]
scaled_predictions = [
    np.exp(p[0].reshape(p[0].shape[0], -1)) * np.exp(p[1]) for p in predictions
]
observed = [np.abs(fold[1]).reshape(fold[1].shape[0], -1) for fold in loci_folds]

log_quantity_correlation = [
    pearsonr(np.log(pred.sum(axis=1) + 1e-3), np.log(obs.sum(axis=1) + 1e-3))[0]
    for pred, obs in zip(scaled_predictions, observed)
]

profile_correlation = [
    pd.DataFrame(pred).corrwith(pd.DataFrame(obs), axis=1)
    for pred, obs in zip(scaled_predictions, observed)
]
