import h5py
import numpy as np
import pandas as pd

oligos = pd.read_csv("media-3_oligos_snps.tsv.gz", sep="\t")
variants = pd.read_csv("media-3_oligos_snps.tsv.gz", sep="\t")["Variant"]

mpra = pd.read_csv("media-4-K562_allelic_mpra.tsv.gz", sep="\t")
mpra["expt"] = np.log2(
    (mpra.mean_RNA_ref_K562 / mpra.mean_Plasmid_ref_K562)
    / (mpra.mean_RNA_alt_K562 / mpra.mean_Plasmid_alt_K562)
)

oligos["pred_p_2114"] = np.log2(
    h5py.File("k562_mpra_snps_2114_ref_procapnet.h5")["quantity"][:]
    / h5py.File("k562_mpra_snps_2114_alt_procapnet.h5")["quantity"][:]
)

oligos["pred_p"] = np.log2(
    h5py.File("../predictions/k562_mpra_snps_ref_procapnet.h5")["quantity"][:]
    / h5py.File("../predictions/k562_mpra_snps_alt_procapnet.h5")["quantity"][:]
)

oligos["pred"] = np.log2(
    h5py.File("../predictions/k562_mpra_snps_ref.h5")["quantity"][:]
    / h5py.File("../predictions/k562_mpra_snps_alt.h5")["quantity"][:]
)

predictions = oligos.merge(mpra, left_on="Variant", right_on="variant")[
    ["Variant", "expt", "pred", "pred_p", "pred_p_2114"]
]
