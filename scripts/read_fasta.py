import gzip
from Bio import SeqIO
import argparse

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

parser = argparse.ArgumentParser()
parser.add_argument("fasta_gz", help="a gzipped (or not) fasta file")
parser.add_argument("out_tsv", help="where to write output tsv")
args = parser.parse_args()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read fasta file to tsv
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def openfile(filename, mode="rt"):
    """
    Handles gzipped files.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


with open(args.out_tsv, "w+") as out:
    out.write("")

with openfile(args.fasta_gz) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = list(record.seq.upper())
        row = [record.id] + sequence
        with open(args.out_tsv, "a") as out:
            out.write("\t".join(row) + "\n")
