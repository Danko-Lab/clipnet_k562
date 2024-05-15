import argparse

import numpy as np
import pandas as pd


def rpm_bedgraph(pl_fn, mn_fn, out_pl_fn, out_mn_fn):
    """Calculates RPM normalization of pl and mn bedgraph files."""
    # read in
    names = ["start", "stop", "counts"]
    pl = pd.read_csv(pl_fn, sep="\t", index_col=0, names=names)
    mn = pd.read_csv(mn_fn, sep="\t", index_col=0, names=names)

    # normalize
    total_reads = pl["counts"].sum() + np.abs(mn["counts"]).sum()
    pl["rpm"] = pl["counts"] * 1e6 / total_reads
    mn["rpm"] = np.abs(mn["counts"]) * 1e6 / total_reads

    # drop counts column
    norm_pl = pl.drop(labels="counts", axis=1)
    norm_mn = mn.drop(labels="counts", axis=1)

    # output
    norm_pl.to_csv(out_pl_fn, sep="\t", header=False)
    norm_mn.to_csv(out_mn_fn, sep="\t", header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pl_fn", help="a (gzipped) bedgraph file (assumes tab separated)"
    )
    parser.add_argument(
        "mn_fn",
        help="a (gzipped) bedgraph file (assumes tab separated). Values can be \
            positive or negative",
    )
    parser.add_argument("out_pl_fn", help="where to write output pl bedgraph file.")
    parser.add_argument("out_mn_fn", help="where to write output mn bedgraph file.")
    args = vars(parser.parse_args())
    rpm_bedgraph(**args)
