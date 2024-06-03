import argparse

import numpy as np
import pandas as pd
import tqdm


def string_to_char_array(seq):
    """
    Converts an ASCII string to a NumPy array of byte-long ASCII codes.
    e.g. "ACGT" becomes [65, 67, 71, 84].
    """
    return np.frombuffer(bytes(seq, "utf8"), dtype=np.int8)


def char_array_to_string(arr):
    """
    Converts a NumPy array of byte-long ASCII codes into an ASCII string.
    e.g. [65, 67, 71, 84] becomes "ACGT".
    """
    assert arr.dtype == np.int8
    return arr.tobytes().decode("ascii")


def one_hot_to_tokens(one_hot):
    """
    Converts an L x D one-hot encoding into an L-vector of integers in the range
    [0, D], where the token D is used when the one-hot encoding is all 0. This
    assumes that the one-hot encoding is well-formed, with at most one 1 in each
    column (and 0s elsewhere).
    """
    tokens = np.tile(one_hot.shape[1], one_hot.shape[0])  # Vector of all D
    seq_inds, dim_inds = np.where(one_hot)
    tokens[seq_inds] = dim_inds
    return tokens


def construct_reporter_seq(enhancer_seq, procapnet=False):
    vector_seqs = {
        "procapnet": "GGCCGAGCGAAGAAGTGGTCCTGCTACTTTGTCCGCCTCCATCCAGTCTATGAGCTGCTGTCGTGATGCTAGAGTAAGAAGTTCGCCAGTGAGTAGTTTCCGAAGAGTTGTGGCCATTGCTACTGGCATCGTGGTATCACGCTCGTCGTTCGGTATGGCTTCGTTCAACTCTGGTTCCCAGCGGTCAAGCCGGGTCACATGATCACCCATATTATGAAGAAATGCAGTCAGCTCCTTAGGGCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCGGTGTTGTCGCTCATGGTAATGGCAGCACTACACAATTCTCTTACCGTCATGCCATCCGTAAGATGCTTTTCCGTGACCGGCGAGTACTCAACCAAGTCGTTTTGTGAGTAGTGTATACGGCGACCAAG",
        "clipnet": "CTGCTCTTGCCCGGCGTCTATACGGGACAACACCGCGCCACATAGCAGTACTTTGAAAGTGCTCATCATCGGGAATCGTTCTTCGGGGCGGAAAGACTCAAGGATCTTGCCGCTATTGAGATCCAGTTCGATATAGCCCACTCTTGCACCCAGTTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCGGGGTGTGCAA",
        "Sfi1_5prime": "AAACAGGCAAGCAAAATGCCGCAAAGAAGGGAATGAGTGCGACACGAAAATGTTGGATGCTCATACTCGTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTACTAGTACGTCTCTCAAGGATAAGTAAGTAATATTAAGGTACGGGAGGTATTGGACAGGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATACGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTCCCCAGTGCAAGTGCAGGTGCCAGAACATTTCTCTGGCCTAAC",
        "reporter_minP_spacer": "",
        "minP_GFP": "tcaatctaaagtatatatgagtaaacttggtctgacagcggccgcaaatgctaaaccactgcagtggttaccagtgcttgatcagtgaggcaccgatctcagcgatctgcctatttcgttcgtccatagtggcctgactccccgtcgtgtagatcactacgattcgtgagggcttaccatcaggccccagcgcagcaatgatgccgcgagagccgcgttcaccggcccccgatttgtcagcaatgaaccagccagcagggagggccgagcgaagaagtggtcctgctactttgtccgcctccatccagtctatgagctgctgtcgtgatgctagagtaagaagttcgccagtgagtagtttccgaagagttgtggccattgctactggcatcgtggtatcacgctcgtcgttcggtatggcttcgttcaactctggttcccagcggtcaagccgggtcacatgatcacccatattatgaagaaatgcagtcagctccttagggcctccgatcgttgtcagaagtaagttggccgcggtgttgtcgctcatggtaatggcagcactacacaattctcttaccgtcatgccatccgtaagatgcttttccgtgaccggcgagtactcaaccaagtcgttttgtgagtagtgtatacggcgaccaagctgctcttgcccggcgtctatacgggacaacaccgcgccacatagcagtactttgaaagtgctcatcatcgggaatcgttcttcggggcggaaagactcaaggatcttgccgctattgagatccagttcgatatagcccactcttgcacccagttgatcttcagcatcttttactttcaccagcgtttcggggtgtgcaaaaacaggcaagcaaaatgccgcaaagaagggaatgagtgcgacacgaaaatgttggatgctcatactcgtcctttttcaatattattgaagcatttatcagggttactagtacgtctctcaag",
    }
    # Join seqs:
    construct = "".join(
        [
            vector_seqs["clipnet"],
            vector_seqs["Sfi1_5prime"][:-5],
            enhancer_seq,
            vector_seqs["reporter_minP_spacer"],
            vector_seqs["minP_GFP"],
        ]
    )
    # Add procapnet sequence if needed
    if procapnet:
        construct = "".join(
            [
                vector_seqs["procapnet"],
                construct,
            ]
        )
        return construct
    # Trim off excess 3' sequence and return encoded sequence
    return construct[:1200]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "oligo_tsv",
        type=str,
        help="Enhancer sequences (oligos) to be inserted into reporter construct. \
            Expects a tsv with column headers 'Variant', 'Allele 1 Oligo', \
            'Allele 2 Oligo', 'chrom', 'pos', 'a1', 'a2'",
    )
    parser.add_argument(
        "reporter_prefix",
        type=str,
        help="Path to reporter construct fasta files to (will write prefix_ref.fa \
            and prefix_alt.fa). If dinuc_shuffles > 0, will write prefix_shuffleN.fa",
    )
    parser.add_argument(
        "--procapnet",
        action="store_true",
        help="Generates reporter constructs for ProCapNet predictions. len = 2114.",
    )
    args = parser.parse_args()
    # Read in oligo sequences
    oligos = pd.read_csv(args.oligo_tsv, sep="\t")
    # Create reporter construct for each oligo sequence and write to file
    with open(f"{args.reporter_prefix}_ref.fa", "w") as ref, open(
        f"{args.reporter_prefix}_alt.fa", "w"
    ) as alt:
        for i, row in tqdm.tqdm(
            list(oligos.iterrows()), desc="Writing reporters", total=oligos.shape[0]
        ):
            ref_construct = construct_reporter_seq(
                row["Allele 1 Oligo"], procapnet=args.procapnet
            )
            alt_construct = construct_reporter_seq(
                row["Allele 2 Oligo"], procapnet=args.procapnet
            )
            name = f"{row['Variant']}_{row['chrom']}_{row['pos']}"
            ref.write(f">{name}_{row['a1']}\n{ref_construct}\n")
            alt.write(f">{name}_{row['a2']}\n{alt_construct}\n")


if __name__ == "__main__":
    main()
