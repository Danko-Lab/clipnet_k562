#!/usr/bin/env python3

import argparse
import pyfastx

# from Bio import SeqIO
import joblib


def chunkify_fasta(fasta_fp, id_fp, output):
    """Extracts one-hot sequences from fasta file by id #."""
    fasta = pyfastx.Fasta(fasta_fp)
    if id_fp is None:
        records = [fasta[i] for i in range(len(fasta))]
    else:
        ids = [int(i.split("_")[0]) - 1 for i in joblib.load(id_fp)]
        records = [fasta[i] for i in ids]
    with open(output, "w") as handle:
        for rec in records:
            handle.write(">%s\n%s\n" % (rec.name, rec.seq))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_fp", help="a fasta file with numerical ids")
    parser.add_argument(
        "--id_fp",
        default=None,
        help="a joblib file w/ list of ids. If not provided, converts entire file",
    )
    parser.add_argument(
        "output", help="where to write output fasta file (uncompressed)"
    )
    args = vars(parser.parse_args())
    chunkify_fasta(**args)
