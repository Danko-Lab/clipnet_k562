### THIS IS BUGGED BECAUSE MP CAUSES ISSUES IN RESAMPLING

import multiprocessing as mp
import random

import numpy as np
import pandas as pd
import tqdm
from scipy.stats import ttest_rel, wilcoxon
from sklearn.metrics import average_precision_score
from sklearn.utils import resample

np.random.seed(47)
random.seed(47)

data = pd.read_csv("/home2/ayh8/clipnet_k562/data/mpra/k562_allelic_mpra_snps.csv.gz")
data = data[
    [
        "fold",
        "variant",
        "emVar_K562",
        "active_K562",
        "log2fc_expt",
        "log2fc_clipnet_ensemble",
        "log2fc_clipnet_holdout",
        "log2fc_procapnet_ensemble",
    ]
]
data = data[np.isfinite(data["log2fc_clipnet_holdout"])]
data.dropna(inplace=True)


def calculate_bootstrap(i):
    emvar, clipnet, procapnet = resample(
        data[data.fold == 0].emVar_K562,
        data[data.fold == 0].log2fc_clipnet_holdout ** 2,
        data[data.fold == 0].log2fc_procapnet_ensemble ** 2,
    )
    return average_precision_score(emvar, clipnet), average_precision_score(
        emvar, procapnet
    )


with mp.Pool(116) as pool:
    results = pd.DataFrame(
        list(tqdm.tqdm(pool.imap(calculate_bootstrap, range(1_000)), total=1_000))
    )

results.columns = ["clipnet", "procapnet"]

wilcoxon(results.clipnet, results.procapnet)
ttest_rel(results.clipnet, results.procapnet)


results.to_csv(
    "/home2/ayh8/clipnet_k562/data/mpra/k562_allelic_mpra_bootstrap_average_precision.csv.gz"
)
