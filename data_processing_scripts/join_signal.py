import argparse
import numpy as np
import pandas as pd


def join_signal(pl, mn, output, in_csv, out_csv, dtype="float16"):
    """joins plus and minus strand signals"""
    if in_csv:
        pl = pd.read_csv(pl, index_col=0, header=None)
        mn = pd.read_csv(mn, index_col=0, header=None)
        joined = pd.concat([pl, mn], axis=1)
    else:
        try:
            pl = np.load(pl)["arr_0"]
            mn = np.load(mn)["arr_0"]
            joined = np.concatenate((pl, mn), axis=-1)
        except ValueError:
            print("Files are not npy format, and is_csv flag not used")
    if out_csv:
        joined.to_csv(output, header=None)
    else:
        if np.array(joined).dtype != dtype:
            joined = joined.astype(dtype)
        np.savez_compressed(output, joined)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pl", type=str, help="a csv or npy file w/ pl signal")
    parser.add_argument("mn", type=str, help="a csv or npy file w/ mn signal")
    parser.add_argument("output", type=str, help="where to write joined output file")
    parser.add_argument(
        "--in_csv",
        help="indicates that input files are csv. default is npy",
        action="store_true",
    )
    parser.add_argument(
        "--out_csv",
        help="indicates that output file is csv. default is npy",
        action="store_true",
    )
    parser.add_argument(
        "--dtype",
        type=str,
        help="what data type the files should be stored as",
        default="float16",
    )
    args = vars(parser.parse_args())

    join_signal(**args)
