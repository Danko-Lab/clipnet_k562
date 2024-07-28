import pandas as pd

# Load sequences

holdout_chroms = list(
    pd.read_csv(
        "/home2/ayh8/clipnet_k562/siraj_mpra/clipnet_data_fold_assignments.csv"
    )[lambda x: x["fold"] == 0]["chrom"]
)

df = pd.read_csv(
    "/home2/ayh8/clipnet_k562/siraj_mpra/media-3_oligos_snps_cleaned.tsv.gz", sep="\t"
)
df["start"] = df["pos"] - 1
df["end"] = df["pos"]
bed = df[["chrom", "start", "end", "Variant"]][
    lambda x: x["chrom"].isin(holdout_chroms)
]
bed.to_csv(
    "/home2/ayh8/clipnet_k562/data/mpra/media-3_oligos_snps_cleaned_holdouts.bed.gz",
    sep="\t",
    index=False,
    header=None,
)

# Then run:

# bedtools slop -l 98304 -r 98303 -i ../data/mpra/media-3_oligos_snps_cleaned_holdouts.bed.gz -g /workdir/ayh8/male.hg19.chrom.sizes | awk 'OFS="\t"{print $1,$2,$3,$4,$3-$2}' | grep 196608 | cut -f1-4 | gzip > ../data/mpra/media-3_oligos_snps_cleaned_enformer_windows_holdouts.bed.gz
# bedtools getfasta -fi /workdir/ayh8/male.hg19.fasta -bed ../data/mpra/media-3_oligos_snps_cleaned_enformer_windows_holdouts.bed.gz -nameOnly | bgzip > ../data/mpra/media-3_oligos_snps_cleaned_enformer_windows_holdouts.fa.gz
