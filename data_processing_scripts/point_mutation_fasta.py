"""
Takes a fasta file of sequences of equal length, introduces a point mutation at a given coordinate.
"""

from Bio import bgzf, Seq, SeqIO, SeqRecord
import pyfastx
import argparse


def point_mutation(seq, position, new_residue):
    """Generate a mutation to new_residue at position in seq"""
    return Seq.Seq(seq[:position] + new_residue + seq[position + 1 :])


def handle(fasta_fp, position, new_residue, out_fp):
    """
    Args
    ---
    fasta_fp    -       filepath of fasta (possibly bgzipped) file. All seqs should be
                        the same length.
    position    -       position where
    new_residue -       new residue
    """
    fna = pyfastx.Fasta(fasta_fp)
    new_records = [
        SeqRecord.SeqRecord(
            point_mutation(rec.seq, position, new_residue), id=rec.name, description=""
        )
        for rec in fna
    ]
    with bgzf.BgzfWriter(out_fp, "wb") as outgz:
        SeqIO.write(sequences=new_records, handle=outgz, format="fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_fp",
        type=str,
        help="a fasta file where all the seqs have the same length.",
    )
    parser.add_argument(
        "position", type=int, help="the position where to introduce a point mutation."
    )
    parser.add_argument(
        "new_residue", type=str, help="the new residue to be introduced at position."
    )
    parser.add_argument(
        "out_fp", type=str, help="where to write output bgzipped fasta file."
    )
    args = vars(parser.parse_args())

    handle(**args)
