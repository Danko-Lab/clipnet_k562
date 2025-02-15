#!/usr/bin/env python3

import argparse
import random

import pandas as pd
import pybedtools
from pybedtools.featurefuncs import midpoint

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser()
parser.add_argument("bed", help="a bed file with genome coordinates")
parser.add_argument(
    "--window",
    help="the size of the windows to get (default = 1000)",
    type=int,
    default=1000,
)
parser.add_argument(
    "--buffer",
    help="min space to window edge on both sides of peak (< window / 2)",
    type=int,
    default=250,
)
parser.add_argument(
    "--reps", help="number of windows to take around each center", type=int, default=1
)
parser.add_argument("--seed", help="a seed for random", default=None)
args = parser.parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform window sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

assert args.window > 0, "Window size needs to be a positive integer."
assert (
    0 <= args.buffer < args.window / 2
), "Buffer needs to be a positive integer less than window / 2."
assert args.reps > 0, "reps needs to be a positive integer."

# set side (for reproducibility)
random.seed(args.seed)

# import data
df = pd.read_csv(args.bed, sep="\t", header=None)


def sample_windows(peaks, window_length, pad, centered=False):
    """Sample windows from self.bed_file"""
    windows = []
    for peak in peaks:
        if centered:
            window_init = midpoint(peak).start
            window_start = window_init - int(window_length / 2)
        else:
            window_init = random.randrange(peak.start, peak.stop + 1)
            window_start = int(window_init - random.randrange(pad, window_length - pad))
        window_stop = window_start + window_length
        windows.append([peak.chrom, window_start, window_stop])
    return pybedtools.BedTool(windows)


# iterate over all rows in df
for i in range(df.shape[0]):
    # perform args.reps samples
    for rep in range(args.reps):
        chr, peak_start, peak_stop = df.iloc[i][:3]
        # pick a random place in the peak to build a window around
        window_init = random.randrange(peak_start, peak_stop + 1)
        # pick a random start point of the window
        window_start = int(
            window_init - random.randrange(args.buffer, args.window - args.buffer)
        )
        window_stop = int(window_start + args.window)
        row = "%s\t%i\t%i" % (chr, window_start, window_stop)
        print(row)
