import os

import pandas as pd

# Download and preprocess MPRA activity data:
os.system(
    "wget https://www.biorxiv.org/content/biorxiv/early/2024/05/06/2024.05.05.592437/DC4/embed/media-4.xlsx"
)
df = pd.read_excel("media-4.xlsx", engine="openpyxl", skiprows=1)
col = ["variant", "active.any", "emVar.any"] + [
    c for c in df.columns if c.endswith("K562")
]
k562 = df[col]
k562.to_csv("media-4-K562_allelic_mpra.tsv.gz", sep="\t", index=False)

# Download and preprocess saturation data:
os.system(
    "wget https://www.biorxiv.org/content/biorxiv/early/2024/05/06/2024.05.05.592437/DC4/embed/media-7.xlsx"
)
saturation = pd.read_excel("media-7.xlsx", engine="openpyxl", skiprows=1)
saturation[saturation["Cell Type"] == "K562"].to_csv(
    "media-7-K562_saturation.tsv.gz", sep="\t", index=False
)

# Download and preprocess oligos:
os.system(
    "wget https://www.biorxiv.org/content/biorxiv/early/2024/05/06/2024.05.05.592437/DC4/embed/media-3.xlsx"
)
oligos = pd.read_excel("media-3.xlsx", engine="openpyxl", skiprows=1)
oligos[["chrom", "pos", "a1", "a2"]] = oligos["Variant"].str.split(":", expand=True)
oligos.to_csv("media-3_oligos.tsv.gz", sep="\t", index=False)
oligos_snps = oligos[(oligos["a1"].str.len() == 1) * (oligos["a2"].str.len() == 1)]
oligos_snps.to_csv("media-3_oligos_snps.tsv.gz", sep="\t", index=False)
