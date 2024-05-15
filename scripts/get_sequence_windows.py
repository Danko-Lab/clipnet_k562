#!/usr/bin/env python3

"""
Takes a directory of fasta files, a directory of bed files, and an output directory.
Uses bedtools getfasta to extract corresponding regions from each fasta file. Assumes
fasta files have a prefix equal matching the first segment of the bed file prefixes.
"""

import argparse
import json
import itertools as it
import multiprocessing as mp
import os


def gsw(exp_dir, exp_prefix, fna_dir, fna_prefix, out_dir):
    bed = os.path.join(exp_dir, "%s_window_uniq.bed" % exp_prefix)
    fna = os.path.join(fna_dir, "%s.fna.gz" % fna_prefix)
    out = os.path.join(out_dir, "%s.fna" % exp_prefix)
    cmd = "pyfastx subseq %s -b %s -o %s" % (fna, bed, out)
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed_dir", type=str, help="input bed directory")
    parser.add_argument("fna_dir", type=str, help="input fasta directory")
    parser.add_argument("out_dir", type=str, help="output directory")
    parser.add_argument(
        "file_prefix_hash",
        type=str,
        help="a json file to convert between experiment and fasta prefix",
    )
    parser.add_argument(
        "--missing",
        type=str,
        default=None,
        help="a list of fasta files that are missing and must be excluded",
    )
    args = parser.parse_args()

    with open(args.file_prefix_hash, "r") as handle:
        d = json.load(handle)
    if args.missing is None:
        missing = []
    else:
        with open(args.missing, "r") as handle:
            missing = handle.read().splitlines()
    exp_prefixes = list(d.keys())
    nonempty_exp_prefixes = [
        prefix for prefix in exp_prefixes if d[prefix] not in missing
    ]
    nonempty_fasta_prefixes = [d[prefix] for prefix in nonempty_exp_prefixes]

    p = mp.Pool(int(mp.cpu_count() / 2))
    p.starmap(
        gsw,
        zip(
            it.repeat(args.bed_dir),
            nonempty_exp_prefixes,
            it.repeat(args.fna_dir),
            nonempty_fasta_prefixes,
            it.repeat(args.out_dir),
        ),
    )
