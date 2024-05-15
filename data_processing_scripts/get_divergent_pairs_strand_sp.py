#!/usr/bin/env python3

import argparse
import pybedtools


def get_divergent_pairs(pl, mn, pl_out, mn_out, distance):
    """Get divergent pairs"""
    pl = pybedtools.BedTool(pl)
    mn = pybedtools.BedTool(mn)

    # sort
    pl_sort = pl.sort()
    mn_sort = mn.sort()

    # get closest with distance
    nb_pl = pl_sort.closest(mn_sort, d=True)
    nb_mn = mn_sort.closest(pl_sort, d=True)

    divergent_pl = []
    divergent_mn = []

    # filter.
    for rec in nb_pl:
        if int(rec[-1]) < distance:
            divergent_pl.append(rec)
    for rec in nb_mn:
        if int(rec[-1]) < distance:
            divergent_mn.append(rec)

    pybedtools.BedTool(divergent_pl).saveas(pl_out)
    pybedtools.BedTool(divergent_mn).saveas(mn_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pl", type=str, help="a bed file containing only pl")
    parser.add_argument("mn", type=str, help="a bed file containing only mn")
    parser.add_argument("pl_out", type=str, help="where to save divergent pl peaks")
    parser.add_argument("mn_out", type=str, help="where to save divergent mn peaks")
    parser.add_argument(
        "--distance", default=200, type=int, help="distance to filter closest pairs on"
    )
    args = vars(parser.parse_args())

    get_divergent_pairs(**args)
