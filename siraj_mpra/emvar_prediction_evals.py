import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.stats import pearsonr, spearmanr, mannwhitneyu, ttest_ind, ks_2samp
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.utils import resample
import multiprocessing as mp

data = pd.read_csv("/home2/ayh8/clipnet_k562/data/mpra/k562_allelic_mpra_snps.csv.gz")

data["pred"] = np.log2(data["ref"] / data["alt"])
data["pred_procapnet"] = np.log2(data["ref_procapnet_ensemble"] / data["alt_procapnet_ensemble"])
for i in range(7):
    data[f"procapnet_fold_{i}"] = np.log2(data[f"ref_procapnet_fold_{i}"] / data[f"alt_procapnet_fold_{i}"])


data = data[np.isfinite(data["pred"])]
data.dropna(inplace=True)

def calculate_bootstrap(i):
    emvar, clipnet, procapnet = resample(
        data[data.fold==0].emVar_K562,
        data[data.fold==0].pred ** 2,
        data[data.fold==0].procapnet_fold_0 ** 2
    )
    return average_precision_score(emvar, clipnet), average_precision_score(emvar, procapnet)

with mp.Pool(50) as pool:
    results = pool.map(calculate_bootstrap, range(int(1e6)))

