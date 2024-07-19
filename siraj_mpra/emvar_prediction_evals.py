import multiprocessing as mp

import numpy as np
import pandas as pd
import tqdm
from sklearn.metrics import average_precision_score
from sklearn.utils import resample

data = pd.read_csv("/home2/ayh8/clipnet_k562/data/mpra/k562_allelic_mpra_snps.csv.gz")

data["pred"] = np.log2(data["ref"] / data["alt"])
data["pred_procapnet"] = np.log2(
    data["ref_procapnet_ensemble"] / data["alt_procapnet_ensemble"]
)
for i in range(7):
    data[f"procapnet_fold_{i}"] = np.log2(
        data[f"ref_procapnet_fold_{i}"] / data[f"alt_procapnet_fold_{i}"]
    )


data = data[np.isfinite(data["pred"])]
data.dropna(inplace=True)


def calculate_bootstrap(i):
    emvar, clipnet, procapnet = resample(
        data[data.fold == 0].emVar_K562,
        data[data.fold == 0].pred ** 2,
        data[data.fold == 0].procapnet_fold_0 ** 2,
    )
    return average_precision_score(emvar, clipnet), average_precision_score(
        emvar, procapnet
    )


with mp.Pool(60) as pool:
    results = list(
        tqdm.tqdm(pool.imap(calculate_bootstrap, range(int(1e6))), total=int(1e6))
    )
