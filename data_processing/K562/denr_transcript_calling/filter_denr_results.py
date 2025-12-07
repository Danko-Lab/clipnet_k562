import numpy as np
import pandas as pd
import pyBigWig
import tqdm


def get_rpk(pybw, chrom, start, end):
    rpb = pybw.stats(chrom, start, end, type="mean")[0]
    return rpb * 1_000 if rpb is not None else 0


df = pd.read_csv(
    "/fs/cbsubscb17/storage/data/hg38/k562/proseq_dreg_datasets/denr_transcript_abundance.bed.gz",
    sep="\t",
    header=None,
)
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

# "/local/workdir/James/PauseEvolution/data/human_K562/K562_LC_1-2_QC_end_all-merge.plus.bw"
with pyBigWig.open("K562_proseq_G156_GC_plus.bw") as bw:
    pl_rpk = [
        get_rpk(bw, row.chrom, row.start, row.end)
        for _, row in tqdm.tqdm(pl_tx.iterrows(), total=pl_tx.shape[0])
    ]

# "/local/workdir/James/PauseEvolution/data/human_K562/K562_LC_1-2_QC_end_all-merge.minus.bw"
with pyBigWig.open("K562_proseq_G156_GC_minus.bw") as bw:
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

df = df[["chrom", "start", "end", "tx_name", "rpb", "strand"]]

df.to_csv(
    "/fs/cbsubscb17/storage/data/hg38/k562/proseq_dreg_datasets/denr_transcript_abundance_rpb.bed.gz",
    sep="\t",
    index=False,
    header=False,
)
df[df["rpb"] >= 1].to_csv(
    "/fs/cbsubscb17/storage/data/hg38/k562/proseq_dreg_datasets/denr_transcript_abundance_rpb_1.bed.gz",
    sep="\t",
    index=False,
    header=False,
)

# Then, run the following script to merge isoforms:
# cd /fs/cbsubscb17/storage/data/hg38/k562/proseq_dreg_datasets/
# gunzip -c denr_transcript_abundance_rpb_1.bed.gz | sort-bed - | bedops -m - | bgzip > denr_greater_than_1rpb_tx.bed.gz
