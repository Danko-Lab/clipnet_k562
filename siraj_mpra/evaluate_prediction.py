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

mpra = pd.read_csv("media-4-K562_allelic_mpra.tsv.gz", sep="\t")
snps = pd.read_csv("media-3_oligos_snps.tsv.gz", sep="\t")
pred = pd.DataFrame({"ref": ref, "alt": alt, "variant": snps["Variant"]})

data = pred.merge(mpra, left_on="variant", right_on="variant")
data["expt"] = np.log2(
    (data["mean_RNA_ref_K562"] / data["mean_Plasmid_ref_K562"])
    / (data["mean_RNA_alt_K562"] / data["mean_Plasmid_alt_K562"])
)
data["pred"] = np.log2(data["ref"] / data["alt"])
data.dropna(inplace=True)

pearsonr(data["expt"], data["pred"])
# PearsonRResult(statistic=0.24444150001598372, pvalue=0.0)