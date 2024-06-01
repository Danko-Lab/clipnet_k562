import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ref = [
    h5py.File(f"/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_fold_{i}.h5")[
        "quantity"
    ][:, 0]
    for i in range(1, 10)
]
alt = [
    h5py.File(f"/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt_fold_{i}.h5")[
        "quantity"
    ][:, 0]
    for i in range(1, 10)
]
ref_p = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_procapnet.h5")[
    "quantity"
]
alt_p = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt_procapnet.h5")[
    "quantity"
]

mpra = pd.read_csv("media-4-K562_allelic_mpra.tsv.gz", sep="\t")
snps = pd.read_csv("media-3_oligos_snps.tsv.gz", sep="\t")

holdouts = pd.read_csv("clipnet_data_fold_assignments.csv")

pred = pd.DataFrame(
    {"ref": ref, "alt": alt, "ref_p": ref_p, "alt_p": alt_p, "variant": snps["Variant"]}
)

data = pred.merge(mpra, left_on="variant", right_on="variant")
data["expt"] = np.log2(
    (data["mean_RNA_ref_K562"] / data["mean_Plasmid_ref_K562"])
    / (data["mean_RNA_alt_K562"] / data["mean_Plasmid_alt_K562"])
)
data["pred"] = np.log2(data["ref"] / data["alt"])
data["pred_p"] = np.log2(data["ref_p"] / data["alt_p"])
data.dropna(inplace=True)

pearsonr(data["expt"], data["pred"])
# PearsonRResult(statistic=0.24444150001598372, pvalue=0.0)
pearsonr(data["expt"], data["pred_p"])
# (0.02508640660560832, 3.072458192675234e-40)
