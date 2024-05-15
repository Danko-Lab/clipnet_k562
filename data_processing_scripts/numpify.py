#!/usr/bin/env python3

import argparse

import joblib
import numpy as np
import pandas as pd


def numpify(in_fn, packbits, output, index_name):
    df = pd.read_csv(in_fn, header=None, index_col=0)
    idx = df.index
    if packbits:
        df = np.packbits(df.astype(int), axis=1)
    np.savez_compressed(output, df)
    joblib.dump(idx, index_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fn", help="a (gzipped) csv file")
    parser.add_argument("output", help="name of output npz file")
    parser.add_argument(
        "--packbits",
        help="apply np.packbits to data? (default = False)",
        action="store_true",
    )
    parser.add_argument("index_name", help="name of index joblib file")
    args = vars(parser.parse_args())
    numpify(**args)
