import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ref = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.h5")["quantity"][
    :, 0
]
alt = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.h5")["quantity"][
    :, 0
]
ref_p = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_procapnet.h5")["quantity"]
alt_p = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt_procapnet.h5")["quantity"]

mpra = pd.read_csv("media-4-K562_allelic_mpra.tsv.gz", sep="\t")
snps = pd.read_csv("media-3_oligos_snps.tsv.gz", sep="\t")
pred = pd.DataFrame({"ref": ref, "alt": alt, "ref_p": ref_p, "alt_p": alt_p, "variant": snps["Variant"]})

data = pred.merge(mpra, left_on="variant", right_on="variant")
data["expt"] = np.log2(
    (data["mean_RNA_ref_K562"] / data["mean_Plasmid_ref_K562"])
    / (data["mean_RNA_alt_K562"] / data["mean_Plasmid_alt_K562"])
)
data["pred"] = np.log2(data["ref"] / data["alt"])
data["pred_p"] = np.log2(data["ref_p"] / data["alt_p"])
data.dropna(inplace=True)

pearsonr(data["expt"], data["pred_p"])
clf = LogisticRegression()
clf.fit(data["pred"].values.reshape(-1, 1), data["emVar_K562"].values)
roc_auc_score(data["emVar_K562"], clf.predict_proba(data["pred"].values.reshape(-1, 1))[:, 1])

# PearsonRResult(statistic=0.24444150001598372, pvalue=0.0)
