#!/usr/bin/env python3

"""
Splits a fasta file into two halves (chr1-10, chr11-22).
"""

from Bio import SeqIO, bgzf
import argparse
import gzip


def openfile(filename, mode="rt"):
    """
    Handles gzipped files.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def main(fasta, first_half, second_half):
    first_match = ["chr%s" % i for i in range(1, 11)]
    second_match = ["chr%s" % i for i in range(11, 23)]
    first_recs = []
    second_recs = []
    with openfile(fasta) as f:
        for rec in SeqIO.parse(f, "fasta"):
            chr = rec.id.split(":")[0].split("_")[-1]
            if chr.replace("rc_", "") in first_match:
                first_recs.append(rec)
            elif chr.replace("rc_", "") in second_match:
                second_recs.append(rec)
            else:
                print("What is this??? %s" % rec.id)
    with bgzf.BgzfWriter(first_half, "wb"):
        SeqIO.write(first_recs, train_out, "fasta")
    with bgzf.BgzfWriter(second_half, "wb") as second_out:
        SeqIO.write(second_recs, second_out, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split the fasta file into individual file with each gene seq"
    )
    parser.add_argument(action="store", dest="fasta", help="Input fasta file")
    parser.add_argument(
        action="store",
        dest="first_half",
        help="Where to output first half (chr1-10) (bgzip)",
    )
    parser.add_argument(
        action="store",
        dest="second_half",
        help="Where to output second half (chr11-22) (bgzip)",
    )
    args = parser.parse_args()
    main(**vars(args))
