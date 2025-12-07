"""
This script calculates the pausing index for a set of genes. The pausing index is defined
as the ratio of the average read density in the promoter region (150 bp upstream and
downstream of the TSS) to the average read density in the gene body (300 bp downstream of
the TSS to 300 bp upstream of the TES). The script takes as input a bed file of gene
coordinates and two bigwig files of read densities (one for the plus strand and one for
the minus strand). The script outputs a bed file with the pausing index for each gene.
"""

import argparse

import numpy as np
import pandas as pd
import pyBigWig

pd.options.mode.copy_on_write = True


def pausing_index(
    bw,
    chrom,
    start,
    stop,
    strand,
    pseudocount=1e-6,
    promoter_boundaries=(150, 150),
    gene_body_boundaries=(300, 300),
):
    """
    Calculate pausing index for a given region.
    IMPORTANT: THIS ASSUMES THAT THE BW FILE IS THE CORRECT STRAND (pausing_index_from_bed
    filters by strand, so this should not be an issue here, but if you use this function
    elsewhere, please check that the supplied bw file is the correct one).
    """
    if strand == "+":
        promoter_boundaries = [
            start - promoter_boundaries[0],
            start + promoter_boundaries[1],
        ]
        gene_body_boundaries = [
            start + gene_body_boundaries[0],
            stop - gene_body_boundaries[1],
        ]
    elif strand == "-":
        promoter_boundaries = [
            stop - promoter_boundaries[0],
            stop + promoter_boundaries[1],
        ]
        gene_body_boundaries = [
            start + gene_body_boundaries[0],
            stop - gene_body_boundaries[1],
        ]
    else:
        raise ValueError("Strand must be either '+' or '-'")
    promoter_read_density = np.abs(
        np.mean(np.nan_to_num(bw.values(chrom, *promoter_boundaries)))
    )
    gene_body_read_density = np.abs(
        np.mean(np.nan_to_num(bw.values(chrom, *gene_body_boundaries)))
    )
    if promoter_read_density is None or gene_body_read_density is None:
        return np.nan
    return promoter_read_density / (gene_body_read_density + pseudocount)


def gene_body_density(
    bw,
    chrom,
    start,
    stop,
    strand,
    gene_body_boundaries=(300, 300),
):
    """
    Calculate pausing index for a given region.
    IMPORTANT: THIS ASSUMES THAT THE BW FILE IS THE CORRECT STRAND (pausing_index_from_bed
    filters by strand, so this should not be an issue here, but if you use this function
    elsewhere, please check that the supplied bw file is the correct one).
    """
    if strand == "+":
        gene_body_boundaries = [
            start + gene_body_boundaries[0],
            stop - gene_body_boundaries[1],
        ]
    elif strand == "-":
        gene_body_boundaries = [
            start + gene_body_boundaries[0],
            stop - gene_body_boundaries[1],
        ]
    else:
        raise ValueError("Strand must be either '+' or '-'")
    gene_body_read_density = np.abs(
        np.mean(np.nan_to_num(bw.values(chrom, *gene_body_boundaries)))
    )
    if gene_body_read_density is None:
        return np.nan
    return gene_body_read_density


def pausing_index_from_bed(
    pl_bw,
    mn_bw,
    bed,
    min_length=2000,
    promoter_boundaries=(150, 150),
    gene_body_boundaries=(300, 300),
    return_gene_body_density=False,
):
    """
    Calculate pausing index for all regions in a bed file. (assumes bed file is a
    pd.DataFrame with at least columns ["chrom", "start", "stop", "name", "strand"]
    """
    bed["gene_length"] = bed["stop"] - bed["start"]
    clean_bed = bed[bed["gene_length"] > min_length]
    pl_bed = clean_bed[clean_bed["strand"] == "+"]
    mn_bed = clean_bed[clean_bed["strand"] == "-"]
    if return_gene_body_density:
        pl_bed["y"] = pl_bed.apply(
            lambda x: gene_body_density(
                pl_bw, x["chrom"], x["start"], x["stop"], x["strand"]
            ),
            axis=1,
        )
        mn_bed["y"] = mn_bed.apply(
            lambda x: gene_body_density(
                mn_bw, x["chrom"], x["start"], x["stop"], x["strand"]
            ),
            axis=1,
        )
        # TPM normalize:
        permillion = 1_000_000 / (pl_bed["y"].sum() + mn_bed["y"].sum())
        pl_bed["y"] = pl_bed["y"] * permillion
        mn_bed["y"] = mn_bed["y"] * permillion
    else:
        pl_bed["y"] = pl_bed.apply(
            lambda x: pausing_index(
                pl_bw, x["chrom"], x["start"], x["stop"], x["strand"]
            ),
            axis=1,
        )
        mn_bed["y"] = mn_bed.apply(
            lambda x: pausing_index(
                mn_bw, x["chrom"], x["start"], x["stop"], x["strand"]
            ),
            axis=1,
        )
    out_bed = pd.concat([pl_bed, mn_bed])
    out_bed["tss"] = np.where(
        out_bed["strand"] == "+", out_bed["start"], out_bed["stop"]
    )
    out_bed["tss+1"] = out_bed["tss"] + 1

    return out_bed[["chrom", "tss", "tss+1", "name", "y", "strand"]]


def main(pl_bw_fp, mn_bw_fp, bed_fp, out_fp, **kwargs):
    """
    Wraps io functions and pausing_index_from_bed.
    """
    bed = pd.read_csv(bed_fp, sep="\t", header=None)
    bed.columns = ["chrom", "start", "stop", "name", "score", "strand"]
    pl_bw = pyBigWig.open(pl_bw_fp)
    mn_bw = pyBigWig.open(mn_bw_fp)
    out_bed = pausing_index_from_bed(pl_bw, mn_bw, bed, **kwargs)
    out_bed.to_csv(out_fp, sep="\t", header=False, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate pausing index for a set of genes."
    )
    parser.add_argument("pl_bw_fp", help="Path to plus strand bigwig file.")
    parser.add_argument("mn_bw_fp", help="Path to minus strand bigwig file.")
    parser.add_argument("bed_fp", help="Path to bed file.")
    parser.add_argument("out_fp", help="Path to output bed file.")
    parser.add_argument(
        "--min_length",
        type=int,
        default=2000,
        help="Minimum gene length to include in output.",
    )
    parser.add_argument(
        "--promoter_boundaries",
        type=int,
        nargs=2,
        default=(150, 150),
        help="Number of bp upstream and downstream of TSS to include in promoter region.",
    )
    parser.add_argument(
        "--gene_body_boundaries",
        type=int,
        nargs=2,
        default=(300, 300),
        help="Number of bp downstream of TSS and upstream of TES to trim off for the gene body.",
    )
    parser.add_argument(
        "--return_gene_body_density",
        action="store_true",
        help="Return gene body read density instead of pausing index.",
    )
    args = parser.parse_args()
    main(
        args.pl_bw_fp,
        args.mn_bw_fp,
        args.bed_fp,
        args.out_fp,
        return_gene_body_density=args.return_gene_body_density,
    )
