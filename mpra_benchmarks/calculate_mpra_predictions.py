import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ref_ensemble = h5py.File(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/clipnet/k562_mpra_snps_ref.h5"
)["quantity"][:, 0]
alt_ensemble = h5py.File(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/clipnet/k562_mpra_snps_alt.h5"
)["quantity"][:, 0]

ref_ref_ensemble = h5py.File(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/clipnet_reference/k562_mpra_snps_ref_ref_model.h5"
)["quantity"][:, 0]

alt_ref_ensemble = h5py.File(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/clipnet_reference/k562_mpra_snps_alt_ref_model.h5"
)["quantity"][:, 0]

ref = [
    h5py.File(
        f"/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/clipnet/k562_mpra_snps_ref_fold_{i}.h5"
    )["quantity"][:, 0]
    for i in range(1, 10)
]
alt = [
    h5py.File(
        f"/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/clipnet/k562_mpra_snps_alt_fold_{i}.h5"
    )["quantity"][:, 0]
    for i in range(1, 10)
]


procapnet_ref = h5py.File(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/procapnet/k562_mpra_snps_2114_ref_procapnet_folds.h5"
)
procapnet_alt = h5py.File(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/procapnet/k562_mpra_snps_2114_alt_procapnet_folds.h5"
)

mpra = pd.read_csv(
    "/home2/ayh8/clipnet_k562/siraj_mpra/media-4-K562_allelic_mpra.tsv.gz",
    sep="\t",
)
snps = pd.read_csv(
    "/home2/ayh8/clipnet_k562/siraj_mpra/media-3_oligos_snps_cleaned.tsv.gz", sep="\t"
)

holdouts = (
    pd.read_csv("/home2/ayh8/clipnet_k562/siraj_mpra/clipnet_data_fold_assignments.csv")
    .set_index("chrom")
    .to_dict()["fold"]
)
holdouts["chrX"] = 0

ref_pred = []
alt_pred = []
folds = []
for i, row in snps.iterrows():
    chrom = row["chrom"]
    if chrom in holdouts:
        fold = holdouts[chrom]
        folds.append(fold)
        if fold == 0:
            ref_pred.append(ref_ensemble[i])
            alt_pred.append(alt_ensemble[i])
        else:
            ref_pred.append(ref[fold - 1][i])
            alt_pred.append(alt[fold - 1][i])

pred = pd.DataFrame(
    {
        "fold": folds,
        "ref_clipnet_ensemble": ref_ensemble,
        "alt_clipnet_ensemble": alt_ensemble,
        "ref_clipnet_reference_ensemble": ref_ref_ensemble,
        "alt_clipnet_reference_ensemble": alt_ref_ensemble,
        "ref_clipnet_holdout": ref_pred,
        "alt_clipnet_holdout": alt_pred,
        "ref_procapnet_ensemble": np.exp(procapnet_ref["ref"][:] - 1).mean(axis=0),
        "alt_procapnet_ensemble": np.exp(procapnet_alt["alt"][:] - 1).mean(axis=0),
        "variant": snps["Variant"],
    }
)


for i in range(7):
    pred[f"ref_procapnet_fold_{i}"] = procapnet_ref["ref"][i]
    pred[f"alt_procapnet_fold_{i}"] = procapnet_alt["alt"][i]

for i in range(1, 10):
    pred[f"ref_clipnet_fold_{i}"] = ref[i - 1]
    pred[f"alt_clipnet_fold_{i}"] = alt[i - 1]

data = pred.merge(mpra, left_on="variant", right_on="variant")
data["log2fc_expt"] = np.log2(
    (data["mean_RNA_ref_K562"] / data["mean_Plasmid_ref_K562"])
    / (data["mean_RNA_alt_K562"] / data["mean_Plasmid_alt_K562"])
)
data["log2fc_clipnet_ensemble"] = np.log2(
    data["ref_clipnet_ensemble"] / data["alt_clipnet_ensemble"]
)
data["log2fc_clipnet_reference_ensemble"] = np.log2(
    data["ref_clipnet_reference_ensemble"] / data["alt_clipnet_reference_ensemble"]
)
data["log2fc_clipnet_holdout"] = np.log2(
    data["ref_clipnet_holdout"] / data["alt_clipnet_holdout"]
)
data["log2fc_procapnet_ensemble"] = np.log2(
    data["ref_procapnet_ensemble"] / data["alt_procapnet_ensemble"]
)

for i in range(7):
    data[f"log2fc_procapnet_fold_{i}"] = np.log2(
        data[f"ref_procapnet_fold_{i}"] / data[f"alt_procapnet_fold_{i}"]
    )

for i in range(1, 10):
    data[f"log2fc_clipnet_fold_{i}"] = np.log2(
        data[f"ref_clipnet_fold_{i}"] / data[f"alt_clipnet_fold_{i}"]
    )

data.to_csv(
    "/fs/cbsubscb17/storage/projects/CLIPNET_transfer/k562_siraj_mpra/k562_allelic_mpra_snps_exp.csv.gz",
    index=False,
)
data = data[np.isfinite(data["log2fc_clipnet_holdout"])]
data.dropna(inplace=True)

pearsons = [
    pearsonr(data[data["fold"] == fold]["expt"], data[data["fold"] == fold]["pred"])[0]
    for fold in range(10)
]

# pearsonr(data[data["fold"] == fold]["expt"], data[data["fold"] == fold]["pred_p"])
pearsonr(data["log2fc_expt"], data["clipnet_pytorch"])
# PearsonRResult(statistic=0.24444150001598372, pvalue=0.0)
pearsonr(data["expt"], data["pred_p"])
# (0.02508640660560832, 3.072458192675234e-40)
