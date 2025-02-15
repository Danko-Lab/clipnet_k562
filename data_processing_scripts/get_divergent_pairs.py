#!/usr/bin/env python3

import argparse
from pybedtools import BedTool


def get_divergent_pairs(pl, mn, distance):
    """Get divergent pairs"""
    pl = BedTool(pl)
    mn = BedTool(mn)

    # sort
    pl_sort = pl.sort()
    mn_sort = mn.sort()

    # get closest with distance
    nb_pl = pl_sort.closest(mn_sort, d=True)
    nb_mn = mn_sort.closest(pl_sort, d=True)

    # filter and print
    for rec in nb_pl:
        if int(rec[-1]) < distance:
            print(str(rec).strip("\n"))

    for rec in nb_mn:
        if int(rec[-1]) < distance:
            print(str(rec).strip("\n"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pl", type=str, help="a bed file containing only pl")
    parser.add_argument("mn", type=str, help="a bed file containing only mn")
    parser.add_argument(
        "--distance", default=200, type=int, help="distance to filter closest pairs on"
    )
    args = vars(parser.parse_args())

    get_divergent_pairs(**args)
