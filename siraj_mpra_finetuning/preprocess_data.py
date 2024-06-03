import pandas as pd
import numpy as np
import pyfastx

chroms = (
    pd.read_csv("clipnet_data_fold_assignments.csv")
    .set_index("chrom")
    .to_dict()["fold"]
)
chroms["chrX"] = 0

mpra = pd.read_csv("../siraj_mpra/media-4-K562_allelic_mpra.tsv.gz", sep="\t")
mpra["chrom"] = mpra["variant"].str.split(":", expand=True)[0]

refs = pyfastx.Fasta("../data/mpra/k562_mpra_snps_ft_ref.fa.gz")
alts = pyfastx.Fasta("../data/mpra/k562_mpra_snps_ft_alt.fa.gz")

snp_names = [x.name.split("_")[0] for x in refs]
ref_seqs = [x.seq for x in refs]
alt_seqs = [x.seq for x in alts]

seqs_df = pd.DataFrame(
    {
        "variant": snp_names,
        "ref_seq": ref_seqs,
        "alt_seq": alt_seqs,
    }
)

data = mpra.merge(seqs_df, on="variant")
data["fold"] = data["chrom"].map(chroms)
data.to_csv("../data/mpra/processed_k562_mpra_data_clipnet_ft.csv.gz", index=False)