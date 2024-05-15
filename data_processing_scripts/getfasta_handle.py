#!/usr/bin/env python3

"""
Takes a directory of fasta files, a directory of bed files, and an output directory.
Uses bedtools getfasta to extract corresponding regions from each fasta file. Assumes
fasta files have a prefix equal matching the first segment of the bed file prefixes.
"""

import argparse
import json
import os


def gsw(fa_dir, bed_fp, out_dir, file_prefix_hash):
    file_prefix = (
        os.path.split(bed_fp)[-1]
        .strip(".bed")
        .strip("_uniq")
        .strip("_intersect_window_test")
    )
    with open(file_prefix_hash, "r") as handle:
        fna_prefix = json.load(handle)[file_prefix]
    fna = os.path.join(fa_dir, "%s.fna" % fna_prefix)
    bed = os.path.join(bed_fp)
    out = os.path.join(out_dir, "%s.fna" % file_prefix)
    cmd = "bedtools getfasta -fi %s -bed %s -fo %s" % (fna, bed, out)
    os.system("echo %s" % cmd)
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fa_dir", type=str, help="fasta directory")
    parser.add_argument("bed_fp", type=str, help="bed file")
    parser.add_argument("out_dir", type=str, help="output directory")
    parser.add_argument(
        "file_prefix_hash",
        type=str,
        help="a json file that converts between procap and fasta prefixes",
    )
    args = parser.parse_args()

    gsw(**vars(args))
