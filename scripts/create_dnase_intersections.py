import argparse
from pybedtools import BedTool

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(description="Extract center positions from a bed file")
parser.add_argument("input", help="filepath to input bed file")
parser.add_argument("output", help="filepath to output txt file")
args = parser.parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import file using BedTool and initialize list of DNase positions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

windows = BedTool(args.input)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loop over all windows and append DNase positions to list
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for i in range(len(windows)):
    w = windows[i]
    start = int(w[1])
    end = int(w[2])
    if int(w[-1]) == end - start:
        entry = (end - start) * ["1"]
    elif int(w[-1]) == 0:
        entry = (end - start) * ["0"]
    else:
        dnase_start = max(start, int(w[-3]))
        dnase_end = min(end, int(w[-2]))
        pre_dnase = dnase_start - start
        dnase = dnase_end - dnase_start
        post_dnase = end - dnase_end
        entry = pre_dnase * ["0"] + dnase * ["1"] + post_dnase * ["0"]
    with open(args.output, "a") as out:
        coords = "%s:%d-%d" % (w[0], start, end)
        out.write(coords + ",".join(entry) + "\n")
