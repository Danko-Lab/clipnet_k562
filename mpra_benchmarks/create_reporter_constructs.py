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


def tokens_to_one_hot(tokens, one_hot_dim):
    """
    Converts an L-vector of integers in the range [0, D] to an L x D one-hot
    encoding. The value `D` must be provided as `one_hot_dim`. A token of D
    means the one-hot encoding is all 0s.
    """
    identity = np.identity(one_hot_dim + 1)[:, :-1]  # Last row is all 0s
    return identity[tokens]


def kshuffle(seq, num_shufs=1, k=2, random_seed=None):
    """
    Creates shuffles of the given sequence, in which dinucleotide frequencies
    are preserved.
    Arguments:
        `seq`: either a string of length L.
        `num_shufs`: the number of shuffles to create, N
        `k`: the length k-mer whose frequencies are to be preserved; defaults
            to k = 2 (i.e. preserve dinucleotide frequencies)
        `rng`: a NumPy RandomState object, to use for performing shuffles
    If `seq` is a string, returns a list of N strings of length L, each one
    being a shuffled version of `seq`.
    """
    # Convert the sequence (string) into a 1D array of numbers (for simplicity)
    if isinstance(seq, str):
        arr = string_to_char_array(seq)
    else:
        raise ValueError("Expected string or one-hot encoded array")

    rng = np.random.RandomState(random_seed)

    if k == 1:
        # Do simple shuffles of `arr`
        all_results = []
        for i in range(num_shufs):
            rng.shuffle(arr)
            all_results.append(char_array_to_string(arr))
        return all_results

    # Tile `arr` from a 1D array to a 2D array of all (k-1)-mers (i.e.
    # "shortmers"), using -1 as a "sentinel" for the last few values
    arr_shortmers = np.empty((len(arr), k - 1), dtype=arr.dtype)
    arr_shortmers[:] = -1
    for i in range(k - 1):
        arr_shortmers[: len(arr) - i, i] = arr[i:]

    # Get the set of all shortmers, and a mapping of which positions start with
    # which shortmers; `tokens` is the mapping, and is an integer representation
    # of the original shortmers (range [0, # unique shortmers - 1])
    shortmers, tokens = np.unique(arr_shortmers, return_inverse=True, axis=0)

    # For each token, get a list of indices of all the tokens that come after it
    shuf_next_inds = []
    for token in range(len(shortmers)):
        # Locations in `arr` where the shortmer exists; some shortmers will have
        # the sentinel, but that's okay
        mask = tokens == token
        inds = np.where(mask)[0]
        shuf_next_inds.append(inds + 1)  # Add 1 to indices for next token

    all_results = []

    for i in range(num_shufs):
        # Shuffle the next indices
        for t in range(len(shortmers)):
            inds = np.arange(len(shuf_next_inds[t]))
            inds[:-1] = rng.permutation(len(inds) - 1)  # Keep last index same
            shuf_next_inds[t] = shuf_next_inds[t][inds]

        counters = [0] * len(shortmers)

        # Build the resulting array
        ind = 0
        result = np.empty_like(tokens)
        result[0] = tokens[ind]
        for j in range(1, len(tokens)):
            t = tokens[ind]
            ind = shuf_next_inds[t][counters[t]]
            counters[t] += 1
            result[j] = tokens[ind]

        shuffled_arr = shortmers[result][:, 0]  # First character of each shortmer
        # (this leaves behind the sentinels)

        all_results.append(char_array_to_string(shuffled_arr))
    return all_results


def construct_reporter_seq(enhancer_seq, procapnet=False):
    vector_seqs = {
        "procapnet_5prime": "GGCCGAGCGAAGAAGTGGTCCTGCTACTTTGTCCGCCTCCATCCAGTCTATGAGCTGCTGTCGTGATGCTAGAGTAAGAAGTTCGCCAGTGAGTAGTTTCCGAAGAGTTGTGGCCATTGCTACTGGCATCGTGGTATCACGCTCGTCGTTCGGTATGGCTTCGTTCAACTCTGGTTCCCAGCGGTCAAGCCGGGTCACATGATCACCCATATTATGAAGAAATGCAGTCAGCTCCTTAGGGCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCGGTGTTGTCGCTCATGGTAATGGCAGCACTACACAATTCTCTTACCGTCATGCCATCCGTAAGATGCTTTTCCGTGACCGGCGAGTACTCAACCAAGTCGTTTTGTGAGTAGTGTATACGGCGACCAAGCTGCTCTTGCCCGGCGTCTATACGGGACAACACCGCGCCACATAGCAGTACTTTGAAAGTGCTCATCATCGGGAATCGTTCTTCGGGGCGGAAAGACTCAAGGATCTTGCCGCTATTGAGATCCAGTTCGATATAGCCCACTCTTGCACCCAGTTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCGGGGTGTGCAA",
        "Sfi1_5prime": "AAACAGGCAAGCAAAATGCCGCAAAGAAGGGAATGAGTGCGACACGAAAATGTTGGATGCTCATACTCGTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTACTAGTACGTCTCTCAAGGATAAGTAAGTAATATTAAGGTACGGGAGGTATTGGACAGGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATACGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTCCCCAGTGCAAGTGCAGGTGCCAGAACATTTCTCTGGCCTAAC",
        "reporter_minP_spacer": "",
        "minP_GFP": "tcaatctaaagtatatatgagtaaacttggtctgacagcggccgcaaatgctaaaccactgcagtggttaccagtgcttgatcagtgaggcaccgatctcagcgatctgcctatttcgttcgtccatagtggcctgactccccgtcgtgtagatcactacgattcgtgagggcttaccatcaggccccagcgcagcaatgatgccgcgagagccgcgttcaccggcccccgatttgtcagcaatgaaccagccagcagggagggccgagcgaagaagtggtcctgctactttgtccgcctccatccagtctatgagctgctgtcgtgatgctagagtaagaagttcgccagtgagtagtttccgaagagttgtggccattgctactggcatcgtggtatcacgctcgtcgttcggtatggcttcgttcaactctggttcccagcggtcaagccgggtcacatgatcacccatattatgaagaaatgcagtcagctccttagggcctccgatcgttgtcagaagtaagttggccgcggtgttgtcgctcatggtaatggcagcactacacaattctcttaccgtcatgccatccgtaagatgcttttccgtgaccggcgagtactcaaccaagtcgttttgtgagtagtgtatacggcgaccaagctgctcttgcccggcgtctatacgggacaacaccgcgccacatagcagtactttgaaagtgctcatcatcgggaatcgttcttcggggcggaaagactcaaggatcttgccgctattgagatccagttcgatatagcccactcttgcacccagttgatcttcagcatcttttactttcaccagcgtttcggggtgtgcaaaaacaggcaagcaaaatgccgcaaagaagggaatgagtgcgacacgaaaatgttggatgctcatactcgtcctttttcaatattattgaagcatttatcagggttactagtacgtctctcaag",
        "procapnet_3prime": "GATAAGTAAGTAATATTAAGGTACGGGAGGTATTGGACAGGCCGCAATAAAATATCTTTATTTTCATTACATCTGTGTGTTGGTTTTTTGTGTGAATCGATAGTACTAACATACGCTCTCCATCAAAACAAAACGAAACAAAACAAACTAGCAAAATAGGCTGTCCCCAGTGCAAGTGCAGGTGCCAGAACATTTCTCTGGCCTAACTGGCCGGTACCTGAGCTCGCTAGCCTCGAGGATATCAAGATCTGGCCTCGGCGGCCAAGCTAGTCGGGGCGGCCGGCCGCTTCGAGCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACAACAACAATTGCATTCATTTTATGTTTCAGGTTAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAAAACCTCTACAAATGTGGTAAAATCGATAAGGATCCGTCGACCGATGCCCTTGAGAGCCTTCAACCCAGTCAGCTCCTTCCGGTGGGCGCGGGGCATGACTATCGTCGCCGCACTTATGACTGTCTTCTTTATCATGCAACTCGTAGGACAGGTGCCGGCAGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGT",
    }
    # Join seqs:
    construct = "".join(
        [
            vector_seqs["Sfi1_5prime"][:-5],
            enhancer_seq,
            vector_seqs["reporter_minP_spacer"],
            vector_seqs["minP_GFP"],
        ]
    )
    # Add procapnet sequence if needed
    if procapnet:
        construct = "".join([vector_seqs["procapnet_5prime"], construct])
        return construct[:2114]
    # Trim off excess 3' sequence and return encoded sequence
    return construct[:1000]


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
        "--dinuc_shuffles",
        type=int,
        default=0,
        help="Number of dinucleotide shuffles to perform on oligo sequences. \
            (if 0, will use original sequence as is).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=47,
        help="Random seed for reproducibility (only used if dinuc_shuffles > 0)",
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
    if args.dinuc_shuffles == 0:
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
    else:
        with open(f"{args.reporter_prefix}_shuffle{args.dinuc_shuffles}.fa", "w") as f:
            for i, row in tqdm.tqdm(
                list(oligos.iterrows()), desc="Writing reporters", total=oligos.shape[0]
            ):
                shuffle = kshuffle(
                    row["Allele 1 Oligo"],
                    num_shufs=args.dinuc_shuffles,
                    random_seed=args.seed,
                )
                for i, s in enumerate(shuffle):
                    construct = construct_reporter_seq(s, procapnet=args.procapnet)
                    f.write(f">{row['Variant']}_shuffle{i}\n{construct}\n")


if __name__ == "__main__":
    main()
