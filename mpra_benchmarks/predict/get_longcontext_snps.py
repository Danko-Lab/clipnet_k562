import pandas as pd

data = pd.read_csv(
    "media-3_oligos_snps_cleaned_holdouts_524288bp_activeK562.bed.gz",
    sep="\t",
    header=None,
)
data = data.iloc[:, 3].str.split(":", expand=True)
data.columns = ["chrom", "pos", "ref", "alt"]
data["pos"] = data["pos"].astype(int) - 1
data.to_csv(
    "media-3_oligos_snps_cleaned_holdouts_524288bp_activeK562.tsv.gz",
    index=False,
    sep="\t",
)
