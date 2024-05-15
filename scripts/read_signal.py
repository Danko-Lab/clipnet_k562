import gzip
import argparse

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser()
parser.add_argument(
    "signal_file", help="a (potentially gzipped) input signal file from bwtools"
)
parser.add_argument(
    "--out", default=None, help="where to write output csv (default: stdout)"
)
args = parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert signal file to csv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def openfile(filename, mode="rt"):
    """Handles gzipped files."""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def parse_signal(signal_str):
    """Helper function that reads a line from a database into a list."""
    signal = signal_str.split(",")
    for i in range(len(signal)):
        if signal[i] in ["NA", "N"]:
            signal[i] = 0
        else:
            signal[i] = abs(float(signal[i]))
    return signal


def get_signal(in_fp, out_fp):
    # iterate over all rows in df
    # extract position and parse signal to readable format
    with openfile(in_fp) as handle:
        if out_fp is not None:
            with open(args.out, "w+") as out:
                out.write("")
        for in_row in handle:
            chrom, start, stop = in_row.split("\t")[0:3]
            signal = parse_signal(in_row.split("\t")[-1].strip("\n"))
            coordinates = f"{chrom}:{str(start)}-{str(stop)}"
            out_row = [coordinates] + signal
            if out_fp is None:
                print(",".join(str(i) for i in out_row))
            else:
                with open(out_fp, "a") as out:
                    out.write(",".join(str(i) for i in out_row) + "\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wrapper
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def main():
    get_signal(args.signal_file, args.out)


if __name__ == "__main__":
    main()
