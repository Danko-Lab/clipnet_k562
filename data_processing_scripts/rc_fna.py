import argparse
from Bio import SeqIO, bgzf
from gzip import open as gzopen


def make_rc_record(record):
    """Returns a new SeqRecord with the reverse complement sequence."""
    return (
        rec.reverse_complement(id="rc_" + rec.id, description="reverse complement")
        for rec in record
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="an input fasta nucleotide file (bgzip)")
    parser.add_argument(
        "output", help="file name for where to write rc fna file (bgzip)."
    )
    args = parser.parse_args()

    rc = make_rc_record(SeqIO.parse(gzopen(args.input, "rt"), "fasta"))
    with bgzf.BgzfWriter(args.output, "wb") as outgz:
        SeqIO.write(rc, outgz, "fasta")
