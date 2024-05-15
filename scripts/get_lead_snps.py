"""
Writes to stdout. Variables named in honor of Billi the cat.
"""

import argparse
import os

import numpy as np
import pandas as pd


def get_lead_snps(input_fp, sig_metric="pvalue", sig_cutoff=1, beta_cutoff=0):
    assert beta_cutoff >= 0, "beta_cutoff must be non-negative."
    assert 0 <= sig_cutoff <= 1, "sig_cutoff must be a probability."
    if input_fp.endswith("csv.gz") or input_fp.endswith("csv"):
        snps = pd.read_csv(input_fp)
    else:
        snps = pd.read_csv(input_fp, sep="\t")
    billi_snps = snps[
        (snps[sig_metric] <= sig_cutoff) & (np.abs(snps["beta"] >= beta_cutoff))
    ]
    billi_names = list(set(list(billi_snps["snps"])))[:100]
    query_snps = ""
    for snp in billi_names:
        query_snps += '"%s"' % snp
        if snp != billi_names[-1]:
            query_snps += ", "
    expr = (
        "'select chrom,chromStart,chromEnd,name,class,convert(alleles using utf8),convert(alleleFreqs using utf8) from snp150 where name in (%s)'"
        % query_snps
    )
    cmd = f"mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -P 3306 -D hg38 --raw --batch -e {expr}"
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fp", help="an input tsv file")
    parser.add_argument(
        "--sig_metric",
        default="pvalue",
        type=str,
        help="what metric to use for signficance",
    )
    parser.add_argument(
        "--sig_cutoff", default=1, type=float, help="a cutoff value for significance"
    )
    parser.add_argument(
        "--beta_cutoff", default=0, type=float, help="a cutoff value for beta"
    )

    get_lead_snps(**vars(parser.parse_args()))
