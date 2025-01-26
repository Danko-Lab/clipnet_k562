import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

ref = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.h5")["quantity"][
    :, 0
]
alt = h5py.File("/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.h5")["quantity"][
    :, 0
]
ref_p = h5py.File(
    "/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_procapnet_folds.h5"
)["ref"][:].mean(axis=0)
alt_p = h5py.File(
    "/home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_procapnet_folds.h5"
)["alt"][:].mean(axis=0)

mpra = pd.read_csv("media-4-K562_allelic_mpra.tsv.gz", sep="\t")
snps = pd.read_csv("media-3_oligos_snps.tsv.gz", sep="\t")
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

pearsonr(data[data.emVar_K562 == 1]["expt"], data[data.emVar_K562 == 1]["pred"])
pearsonr(data[data.emVar_K562 == 1]["expt"], data[data.emVar_K562 == 1]["pred_p"])
