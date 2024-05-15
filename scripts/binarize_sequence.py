#!/usr/bin/env python

import numpy as np
import pandas as pd
import argparse
from joblib import dump


def process(sequence_df):
    X = []
    for i in range(sequence_df.shape[0]):
        sequence = sequence_df.iloc[i,]
        sequence_na = sequence.replace({"N": np.nan})  # Mark Ns as NaNs
        X.append(pd.get_dummies(sequence_na.transpose()))
    return X


def main(in_fn, out_fn):
    sequence_df = pd.read_csv(in_fn, sep="\t", header=None, index_col=0)
    output_df = process(sequence_df)
    dump(output_df, out_fn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("in_fn", help="a tsv file with window sequence")
    parser.add_argument("out_fn", help="where to write output tsv")
    args = vars(parser.parse_args())

    main(**args)
