import numpy as np
import pandas as pd
import pyBigWig
import tqdm


def get_rpk(pybw, chrom, start, end):
    rpb = pybw.stats(chrom, start, end, type="mean")[0]
    return rpb * 1_000 if rpb is not None else 0


df = pd.read_csv("denr_transcript_abundance.bed.gz", sep="\t", header=None)
df.columns = [
    "chrom",
    "start",
    "end",
    "strand",
    "tx_name",
    "tx_abundance",
    "tx_id",
    "not_sure",
]
df["start"] = df.start - 1

pl_tx = df[df.strand == "+"]
mn_tx = df[df.strand == "-"]

with pyBigWig.open(
    "/local/workdir/James/PauseEvolution/data/human_K562/K562_LC_1-2_QC_end_all-merge.plus.bw"
) as bw:
    pl_rpk = [
        get_rpk(bw, row.chrom, row.start, row.end)
        for _, row in tqdm.tqdm(pl_tx.iterrows(), total=pl_tx.shape[0])
    ]

with pyBigWig.open(
    "/local/workdir/James/PauseEvolution/data/human_K562/K562_LC_1-2_QC_end_all-merge.minus.bw"
) as bw:
    mn_rpk = [
        get_rpk(bw, row.chrom, row.start, row.end)
        for _, row in tqdm.tqdm(mn_tx.iterrows(), total=mn_tx.shape[0])
    ]

pl_tx["rpk"] = pl_rpk
mn_tx["rpk"] = np.abs(mn_rpk)

df = pd.concat([pl_tx, mn_tx])
df["tx_abundance"] > 0

per_million = df.rpk.sum() / 1_000_000_000

df["rpb"] = df.tx_abundance * per_million

df.to_csv("denr_transcript_abundance_rpb.bed.gz", sep="\t", index=False, header=False)
df[df["rpb"] >= 1].to_csv(
    "denr_transcript_abundance_rpb_1.bed.gz", sep="\t", index=False, header=False
)

# Then, run the following script to merge isoforms:
# gunzip -c denr_transcript_abundance_rpb_1.bed.gz | sort-bed - | bedtools merge | bgzip > denr_greater_than_1rpb_tx.bed.gz
