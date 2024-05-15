import argparse

import pandas as pd

autosomes = ["chr%d" % i for i in range(1, 23)]


def clean(input_fp):
    df = pd.read_csv(input_fp, sep="\t")
    autosome_only = df[df["chrom"].isin(autosomes)].reset_index(drop=True)
    cleaned = (
        autosome_only.dropna()
        .drop_duplicates("chromStart", keep=False)
        .reset_index(drop=True)
    )
    binary_only = cleaned[
        [
            cleaned["class"][i] == "single"
            and len(cleaned["convert(alleles using utf8)"][i]) == 4
            for i in range(cleaned.shape[0])
        ]
    ].reset_index(drop=True)
    return binary_only


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fp", type=str, help="an input tsv with a chrom column")
    parser.add_argument("output_fp", type=str, help="where to write output tsv")
    args = parser.parse_args()
    output = clean(args.input_fp)
    output.to_csv(args.output_fp, sep="\t", index=False)
