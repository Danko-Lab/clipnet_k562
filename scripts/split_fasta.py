#!/usr/bin/env python3

from Bio import SeqIO, bgzf
import argparse
import gzip

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser(
    description="Split the fasta file into individual file with each gene seq"
)
parser.add_argument(action="store", dest="input", help="Input fasta file")
parser.add_argument(
    action="store",
    dest="train_out",
    help="Where to output train fasta file (chr1-19) (bgzip)",
)
parser.add_argument(
    action="store",
    dest="val_out",
    help="Where to output validation fasta file (chr20) (bgzip)",
)
parser.add_argument(
    action="store",
    dest="test_out",
    help="Where to output test fasta file (chr21,22) (bgzip)",
)
args = parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Split fasta into training and testing sets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def openfile(filename, mode="rt"):
    """
    Handles gzipped files.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


train_match = ["chr%s" % i for i in range(1, 20)]
val_match = ["chr20"]
test_match = ["chr%s" % i for i in range(21, 23)]

train_recs = []
val_recs = []
test_recs = []

# with openfile(args.input) as f:
#    for rec in SeqIO.parse(f, 'fasta'):
#        chr = rec.id.split('_')[1]
#        if chr.id.split(':')[0].replace('rc_', '') in train_match:
#            train_recs.append(rec)
#        elif chr.id.split(':')[0].replace('rc_', '') in val_match:
#            val_recs.append(rec)
#        elif chr.id.split(':')[0].replace('rc_', '') in test_match:
#            test_recs.append(rec)

with openfile(args.input) as f:
    for rec in SeqIO.parse(f, "fasta"):
        chr = rec.id.split(":")[0].split("_")[-1]
        if chr.replace("rc_", "") in train_match:
            train_recs.append(rec)
        elif chr.replace("rc_", "") in val_match:
            val_recs.append(rec)
        elif chr.replace("rc_", "") in test_match:
            test_recs.append(rec)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write to output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with bgzf.BgzfWriter(args.train_out, "wb") as train_out:
    SeqIO.write(train_recs, train_out, "fasta")

with bgzf.BgzfWriter(args.val_out, "wb") as val_out:
    SeqIO.write(val_recs, val_out, "fasta")

with bgzf.BgzfWriter(args.test_out, "wb") as test_out:
    SeqIO.write(test_recs, test_out, "fasta")
