#!/usr/bin/env python3

import multiprocessing as mp
import os

with open(
    "/fs/cbsudanko/storage/projects/CLIPNET/clipnet_scripts/preprocessing_scripts/procap_file_prefixes.txt"
) as f:
    file_prefixes = f.read().splitlines()

WORKDIR = "/workdir/ayh8/gm/"
DANKO_0001 = "/home/danko_0001/projects/ayh8/"


def gsw(file_prefix):
    individual_id = file_prefix.split("_")[0]
    fna = os.path.join(DANKO_0001, "1000genomes_yrb/consensus/%s.fna" % individual_id)
    bed = os.path.join(PATH, "nonqtl/nonvariable_expression_random_windows_uniq.bed")
    out = os.path.join(PATH, "nonqtl/sequence/%s.fna" % file_prefix)
    cmd = "bedtools getfasta -fi %s -bed %s -fo %s" % (fna, bed, out)
    os.system(cmd)


if __name__ == "__main__":
    p = mp.Pool(int(mp.cpu_count() / 2))
    p.map(gsw, file_prefixes)
