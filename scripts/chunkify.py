#!/usr/bin/env python3

import argparse
import joblib
import numpy as np
import pandas as pd


def chunkify(in_fn, output_prefix, packbits=False, chunksize=2**18, dtype="float16"):
    """Breaks a csv file into several chunks."""
    assert dtype in [
        "float16",
        "float32",
        "float64",
    ], "dtype must be float16, float32, or float64"
    n_chunk = 0
    for chunk in pd.read_csv(in_fn, header=None, index_col=0, chunksize=chunksize):
        chunkname = "%s_%d" % (output_prefix, n_chunk)
        idx = chunk.index
        if packbits:
            chunk = np.packbits(chunk.astype(int), axis=1)
        else:
            chunk = chunk.astype(dtype)
        print("Writing chunk %d data to %s.npz ... " % (n_chunk, chunkname))
        np.savez_compressed("%s.npz" % chunkname, chunk)
        joblib.dump(idx, "%s_idx.joblib.gz" % chunkname)
        n_chunk += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fn", type=str, help="a (gzipped) csv file")
    parser.add_argument("output_prefix", type=str, help="prefix for output files")
    parser.add_argument(
        "--packbits",
        help="apply np.packbits to data? (default = False)",
        action="store_true",
    )
    parser.add_argument(
        "--chunksize",
        help="size of chunks to break files into",
        type=int,
        default=2**18,
    )
    parser.add_argument(
        "--dtype",
        help="what dtype to save the array as.",
        type=str,
        default="float16",
    )
    args = vars(parser.parse_args())
    chunkify(**args)
