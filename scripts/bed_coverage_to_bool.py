#!/usr/bin/env python3

"""
Takes a bed file of windows (fixed length, longer regions) and a bed file of peaks
(variable length, representing locations of some signal, shorter regions). Outputs a
dataframe with boolean values indicating whether each position in each window intersects
with a peak.
"""

import argparse
import itertools
import multiprocessing as mp
import pandas as pd
import pybedtools
import numpy as np


def get_coverage(windows, peaks, i):
    window = tuple(windows[i][0:3])
    name = "%s:%s-%s" % window
    coverage = pybedtools.BedTool("%s %s %s" % window, from_string=True).coverage(
        peaks, d=True
    )
    return name, coverage.to_dataframe()["score"]


def bed_coverage_to_bool_mp(peaks_bed, windows_bed, out_csv):
    windows = pybedtools.BedTool(windows_bed)
    peaks = pybedtools.BedTool(peaks_bed)
    indices = range(len(windows))
    with mp.Pool(processes=max(16, int(mp.cpu_count() / 2))) as pool:
        coverage = pool.starmap(
            get_coverage,
            zip(itertools.repeat(windows), itertools.repeat(peaks), indices),
        )
    df = pd.DataFrame(coverage).transpose()
    df.to_csv(out_csv)


def bed_coverage_to_bool(peaks_bed, windows_bed, out_csv):
    windows = pybedtools.BedTool(windows_bed)
    peaks = pybedtools.BedTool(peaks_bed).sort().merge()
    boolean = {}
    for i in range(len(windows)):
        window = tuple(windows[i][0:3])
        name = "%s:%s-%s" % window
        coverage = pybedtools.BedTool("%s %s %s" % window, from_string=True).coverage(
            peaks, d=True
        )
        boolean[name] = np.clip(coverage.to_dataframe()["score"], 0, 1)
    df = pd.DataFrame(boolean).transpose()
    df.to_csv(out_csv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "windows_bed", type=str, help="where to load windows from (longer than peaks)"
    )
    parser.add_argument(
        "peaks_bed", type=str, help="where to load peaks from (shorter than windows)"
    )
    parser.add_argument("out_csv", type=str, help="name of csv to output to")
    args = parser.parse_args()

    bed_coverage_to_bool(**vars(args))
