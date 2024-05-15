"""
Orient fasta and procap files so that the max TSS is on the plus strand.
"""

import argparse
import os
import numpy as np
import pandas as pd
from Bio import bgzf, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyfastx


def load_track(fp, unpackbits=False):
    if os.path.splitext(fp)[-1] == ".npz":
        arr = np.load(fp)["arr_0"]
        if unpackbits:
            return np.unpackbits(arr, axis=1)
        else:
            return arr
    else:
        return np.array(pd.read_csv(fp, index_col=0, header=None))


def calc_tss_strand(fasta_fp, procap_fp):
    fasta = pyfastx.Fasta(fasta_fp)
    procap = load_track(procap_fp)
    tss_pl = np.argmax(procap[:, 250:750], axis=1)
    tss_mn = np.argmax(procap[:, 1250:1750], axis=1)
    flip = [False if tss_pl[i] >= tss_mn[i] else True for i in range(tss_pl.shape[0])]
    out_fasta = [
        SeqRecord(Seq(fasta[i].antisense), id=str(fasta[i].id), name=fasta[i].name)
        if flip[i]
        else SeqRecord(Seq(fasta[i].seq), id=str(fasta[i].id), name=fasta[i].name)
        for i in range(len(flip))
    ]
    out_procap = np.array(
        [procap[i, ::-1] if flip[i] else procap[i, :] for i in range(len(flip))]
    )
    return out_fasta, out_procap


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_fp", type=str, help="a file path to the fasta file.")
    parser.add_argument(
        "procap_fp",
        type=str,
        help="a file path to the procap file (used to calculate tss positions).",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
        help="where to write output. Will write in same directory if left blank.",
    )
    args = parser.parse_args()

    if args.out_dir is None:
        out_dir = os.path.split(args.fasta_fp)[0]
    else:
        out_dir = args.out_dir
    fasta_fn = os.path.split(args.fasta_fp)[-1]
    out_fasta_fp = os.path.join(
        out_dir, "%s_oriented.fna.gz" % fasta_fn.partition(".")[0]
    )
    procap_fn = os.path.split(args.procap_fp)[-1]
    out_procap_fp = os.path.join(
        out_dir, "%s_oriented.csv.gz" % procap_fn.partition(".")[0]
    )

    output = calc_tss_strand(args.fasta_fp, args.procap_fp)
    with bgzf.BgzfWriter(out_fasta_fp, "wb") as outgz:
        SeqIO.write(sequences=output[0], handle=outgz, format="fasta")
    pd.DataFrame(output[1]).to_csv(out_procap_fp, header=None)
