#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2019-11-06
#
# PINTS: Peak Identifier for Nascent Transcripts Sequencing
# Copyright (C) 2019 Li Yao at the Yu Lab
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
import logging
import os
import sys

import numpy as np
import pandas as pd
import pysam

NP_DT_RANGE = (127, 32767, 2147483647, 9223372036854775807)
NP_DT_NAME = (np.int8, np.int16, np.int32, np.int64)


def log_assert(bool_, message, logger):
    """
    Assert function, result will be wrote to logger

    Parameters
    ----------
    bool_ : bool expression
        Condition which should be true
    message : str
        This message will show up when bool_ is False
    logger : logger
        Specific logger

    Returns
    -------

    """
    try:
        assert bool_, message
    except AssertionError as err:
        logger.error("%s" % err)
        sys.exit(1)


def _get_coverage_bw(bw_file, chromosome_startswith, output_dir, output_prefix, **kwargs):
    """
    Get coverage information from bigwig files

    Parameters
    ----------
    bw_file : str
        Path to bigwig file
    chromosome_startswith : str
         Filter out reads whose reference don't start with chromosome_startswith
    output_dir : str
        Path to write outputs
    output_prefix : str
        Prefix of outputs
    kwargs :

    Returns
    -------

    """
    logger = logging.getLogger("PINTS - IO engine")
    try:
        import pyBigWig
    except ImportError as e:
        missing_package = str(e).replace("No module named '", "").replace("'", "")
        logger.error("Please install %s first!" % missing_package)
        sys.exit(-1)
    bw = pyBigWig.open(bw_file)
    log_assert(bw.isBigWig(), "BigWig file %s is not valid." % bw_file, logger)
    chromosome_coverage = dict()
    chromosome_pre_accessible_dict = dict()
    chromosome_reads_dict = dict()

    data_type = np.int32
    bw_max = max(abs(bw.header()["maxVal"]), abs(bw.header()["minVal"]))
    for k, v in enumerate(NP_DT_RANGE):
        if v > bw_max:
            data_type = NP_DT_NAME[k]
            break

    chrom_size_list = []
    # max_reference_length = 0
    genome_size = 0

    # result = {}
    for chromosome, csize in bw.chroms().items():
        if chromosome.startswith(chromosome_startswith):
            _chromosome_cov = np.nan_to_num(bw.values(chromosome, 0, csize))
            _chromosome_cov[_chromosome_cov < 0] *= -1
            _chromosome_cov = _chromosome_cov.astype(data_type, copy=False)
            fn = os.path.join(output_dir, "%s_%s" % (output_prefix, chromosome))
            np.save(fn, _chromosome_cov)
            chromosome_coverage[chromosome] = fn + ".npy"
            chrom_size_list.append((chromosome, csize))
            chromosome_pre_accessible_dict[chromosome] = []
            chromosome_reads_dict[chromosome] = 0
            genome_size += csize
    """
    if "output_chrom_size" in kwargs.keys() and kwargs["output_chrom_size"]:
        csize_fn = os.path.join(output_dir, output_prefix + ".csize")
        with open(csize_fn, "w") as f:
            pd.DataFrame(chrom_size_list, columns=["Chromosome", "Size"]).to_csv(f,
                                                                                 sep="\t",
                                                                                 index=False,
                                                                                 header=False)
    """
    bw.close()
    return chromosome_coverage


def get_coverage_bw(bw_pl, bw_mn, chromosome_startswith, output_dir, output_prefix, **kwargs):
    """
    Get 5' coverage

    Parameters
    ----------
    bw_pl : str
        Full path to input bigwig file, plus strand
    bw_mn : str
        Full path to input bigwig file, minus strand
    chromosome_startswith : str
        Filter out reads whose reference don't start with chromosome_startswith
    output_dir : str
        Path to write outputs
    output_prefix : str
        Prefix of outputs
    **kwargs


    Returns
    -------
    pl : dict
        Dictionary of per base coverage per chromosome (positive strand)
    mn : dict
        Dictionary of per base coverage per chromosome (negative strand)
    """
    logger = logging.getLogger("IO engine")
    pl = _get_coverage_bw(bw_pl, chromosome_startswith, output_dir, output_prefix+"_pl", **kwargs)
    mn = _get_coverage_bw(bw_mn, chromosome_startswith, output_dir, output_prefix+"_mn", **kwargs)
    log_assert(pl.keys() == mn.keys(), "bw_pl and bw_mn should have the same chromosomes", logger)
    return pl, mn


def _bam_se_5p_ss(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    """

    Parameters
    ----------
    bam_obj
    pl_cov
    mn_cov
    chromosome_startswith
    mapq
    kwargs

    Returns
    -------

    """
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1

        if read_strand == "+":
            pl_cov[read.reference_name][pos_5prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_5prime] += 1


def _bam_se_5p_rc(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.reference_name not in chroms:
            continue
        if read.mapq < mapq:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1
        read_strand = "-" if read_strand == "+" else "+"

        if read_strand == "+":
            pl_cov[read.reference_name][pos_5prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_5prime] += 1


def _bam_se_3p_ss(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_3prime = read.reference_end - 1 if read_strand == "+" else read.reference_start

        if read_strand == "+":
            pl_cov[read.reference_name][pos_3prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_3prime] += 1


def _bam_se_3p_rc(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_3prime = read.reference_end - 1 if read_strand == "+" else read.reference_start
        read_strand = "-" if read_strand == "+" else "+"

        if read_strand == "+":
            pl_cov[read.reference_name][pos_3prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_3prime] += 1


def _bam_pe_r1_5p_ss(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read2 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1

        if read_strand == "+":
            pl_cov[read.reference_name][pos_5prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_5prime] += 1


def _bam_pe_r1_5p_rc(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read2 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1
        read_strand = "-" if read_strand == "+" else "+"

        if read_strand == "+":
            pl_cov[read.reference_name][pos_5prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_5prime] += 1


def _bam_pe_r2_5p_ss(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read1 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1

        if read_strand == "+":
            pl_cov[read.reference_name][pos_5prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_5prime] += 1


def _bam_pe_r2_5p_rc(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read1 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1
        read_strand = "-" if read_strand == "+" else "+"

        if read_strand == "+":
            pl_cov[read.reference_name][pos_5prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_5prime] += 1


def _bam_pe_r1_3p_ss(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read2 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_3prime = read.reference_end - 1 if read_strand == "+" else read.reference_start

        if read_strand == "+":
            pl_cov[read.reference_name][pos_3prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_3prime] += 1


def _bam_pe_r1_3p_rc(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read2 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_3prime = read.reference_end - 1 if read_strand == "+" else read.reference_start
        read_strand = "-" if read_strand == "+" else "+"

        if read_strand == "+":
            pl_cov[read.reference_name][pos_3prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_3prime] += 1


def _bam_pe_r2_3p_ss(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read1 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_3prime = read.reference_end - 1 if read_strand == "+" else read.reference_start

        if read_strand == "+":
            pl_cov[read.reference_name][pos_3prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_3prime] += 1


def _bam_pe_r2_3p_rc(bam_obj, pl_cov, mn_cov, mapq, **kwargs):
    c = list(pl_cov.keys())
    c.extend(list(mn_cov.keys()))
    if c is not None:
        chroms = set(c)
    else:
        chroms = set()
    for read in bam_obj:
        if read.is_read1 or read.mapq < mapq:
            continue
        if read.reference_name not in chroms:
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_3prime = read.reference_end - 1 if read_strand == "+" else read.reference_start
        read_strand = "-" if read_strand == "+" else "+"

        if read_strand == "+":
            pl_cov[read.reference_name][pos_3prime] += 1
        elif read_strand == "-":
            mn_cov[read.reference_name][pos_3prime] += 1


def get_read_signal(input_bam, loc_prime, chromosome_startswith, output_dir, output_prefix, filters=[],
                    reverse_complement=False, **kwargs):
    """

    Parameters
    ----------
    input_bam : input_bam : str
        Full path to input bam file
    loc_prime : str, format: Read_End

    chromosome_startswith : str or None
        Filter out reads whose reference don't start with chromosome_startswith
    output_dir : str
        Path to write outputs
    output_prefix : str
        Prefix of outputs
    filters : list or tuple
        List of keywords to filter chromosomes
    reverse_complement :
    **kwargs :
        Keyword arguments, mapq_threshold (int, default 30) and output_chrom_size (bool, default True) is effective
    Returns
    -------
    chromosome_coverage_pl : dict
        Dictionary of per base coverage per chromosome (positive strand)
    chromosome_coverage_mn : dict
        Dictionary of per base coverage per chromosome (negative strand)
    """
    logger = logging.getLogger("IO engine")
    supported_protocols = {
        "PROcap": "R_5_f",
        "GROcap": "R_5_f",
        "PROseq": "R_5_r",
        "GROseq": "R_5_f",
        "CoPRO": "R2_5_f",
        "CAGE": "R_5_f",
    }
    if loc_prime in supported_protocols:
        read_num, read_end, fw_rv = supported_protocols[loc_prime].split("_")
        loc_prime = "%s_%s" % (read_num, read_end)
        reverse_complement = False if fw_rv == "f" else True
    log_assert(loc_prime in ("R_5", "R_3", "R1_5", "R1_3", "R2_5", "R2_3"),
               "library_type must be R1_5, R1_3, R2_5 or R2_3", logger)
    library_layout, interested_end = loc_prime.split("_")
    mapq_threshold = kwargs["mapq_threshold"] if "mapq_threshold" in kwargs.keys() else 30
    chromosome_coverage_pl = dict()
    chromosome_coverage_mn = dict()

    result_pl = {}
    result_mn = {}

    try:
        bam = pysam.AlignmentFile(input_bam, "rb")
    except Exception as e:
        logger.error(e)
        sys.exit(-1)
    genome_size = 0

    for chromosome in bam.header["SQ"]:
        if "SN" in chromosome.keys() and "LN" in chromosome.keys() and \
                chromosome["SN"].startswith(chromosome_startswith):
            ignore_flag = 0
            for keyword in filters:
                if chromosome["SN"].find(keyword) != -1:
                    ignore_flag = 1
            if not ignore_flag:
                chromosome_coverage_pl[chromosome["SN"]] = np.zeros(chromosome["LN"], dtype=np.uint32)
                chromosome_coverage_mn[chromosome["SN"]] = np.zeros(chromosome["LN"], dtype=np.uint32)
                genome_size += chromosome["LN"]

    load_cache_flag = 1
    for chrom in chromosome_coverage_pl:
        fn = os.path.join(output_dir, output_prefix + "_pl_%s" % chrom) + ".npy"
        if not os.path.exists(fn):
            load_cache_flag = 0
            break
        result_pl[chrom] = fn
    for chrom in chromosome_coverage_mn:
        fn = os.path.join(output_dir, output_prefix + "_mn_%s" % chrom) + ".npy"
        if not os.path.exists(fn):
            load_cache_flag = 0
            break
        result_mn[chrom] = fn
    if load_cache_flag:
        logger.info("Loading cache from previous result, if you want to re-parse the bam file, "
                    "please delete all files end with npy first.")
    else:
        bam_parser = None
        if library_layout == "R":
            if interested_end == "5":
                if reverse_complement:
                    bam_parser = _bam_se_5p_rc
                else:
                    bam_parser = _bam_se_5p_ss
            else:
                if reverse_complement:
                    bam_parser = _bam_se_3p_rc
                else:
                    bam_parser = _bam_se_3p_ss
        elif library_layout == "R1":
            if interested_end == "5":
                if reverse_complement:
                    bam_parser = _bam_pe_r1_5p_rc
                else:
                    bam_parser = _bam_pe_r1_5p_ss
            else:
                if reverse_complement:
                    bam_parser = _bam_pe_r1_3p_rc
                else:
                    bam_parser = _bam_pe_r1_3p_ss
        elif library_layout == "R2":
            if interested_end == "5":
                if reverse_complement:
                    bam_parser = _bam_pe_r2_5p_rc
                else:
                    bam_parser = _bam_pe_r2_5p_ss
            else:
                if reverse_complement:
                    bam_parser = _bam_pe_r2_3p_rc
                else:
                    bam_parser = _bam_pe_r2_3p_ss

        log_assert(bam_parser is not None, "Cannot initiate a parser for this experiment", logger)

        bam_parser(bam_obj=bam, pl_cov=chromosome_coverage_pl, mn_cov=chromosome_coverage_mn, mapq=mapq_threshold)

        # convert data type if necessary
        for chrom in chromosome_coverage_pl.keys():
            pl_chrom_max = np.max(chromosome_coverage_pl[chrom])
            mn_chrom_max = np.max(chromosome_coverage_mn[chrom])
            for k, v in enumerate(NP_DT_RANGE):
                if v > pl_chrom_max:
                    chromosome_coverage_pl[chrom] = chromosome_coverage_pl[chrom].astype(NP_DT_NAME[k], copy=False)
                    fn = os.path.join(output_dir, output_prefix + "_pl_%s" % chrom)
                    np.save(fn, chromosome_coverage_pl[chrom])
                    result_pl[chrom] = fn+".npy"
                    break
            for k, v in enumerate(NP_DT_RANGE):
                if v > mn_chrom_max:
                    chromosome_coverage_mn[chrom] = chromosome_coverage_mn[chrom].astype(NP_DT_NAME[k], copy=False)
                    fn = os.path.join(output_dir, output_prefix + "_mn_%s" % chrom)
                    np.save(fn, chromosome_coverage_mn[chrom])
                    result_mn[chrom] = fn + ".npy"
                    break
    return result_pl, result_mn


def get_coverage(input_bam, library_type, chromosome_startswith, output_dir, output_prefix, **kwargs):
    """
    Get 5' coverage

    Parameters
    ----------
    input_bam : str
        Full path to input bam file
    library_type : str
        Library type, available options: se, pe_fr, pe_rf
    chromosome_startswith : str
        Filter out reads whose reference don't start with chromosome_startswith
    output_dir : str
        Path to write outputs
    output_prefix : str
        Prefix of outputs
    **kwargs
        Keyword arguments, mapq_threshold (int, default 30) and output_chrom_size (bool, default True) is effective

    Returns
    -------
    chromosome_coverage_pl : dict
        Dictionary of per base coverage per chromosome (positive strand)
    chromosome_coverage_mn : dict
        Dictionary of per base coverage per chromosome (negative strand)
    """
    logger = logging.getLogger("IO engine")
    log_assert(library_type in ("se", "pe_fr", "pe_rf"), "library_type must be se, pe_fr or pe_rf", logger)
    mapq_threshold = kwargs["mapq_threshold"] if "mapq_threshold" in kwargs.keys() else 30
    chromosome_coverage_pl = dict()
    chromosome_coverage_mn = dict()
    chromosome_pre_accessible_dict = dict()
    chromosome_reads_dict = dict()
    # chromosome_lambda_dict = dict()
    bam = pysam.AlignmentFile(input_bam, "rb")

    chrom_size_list = []
    # max_reference_length = 0
    genome_size = 0

    for chromosome in bam.header["SQ"]:
        if "SN" in chromosome.keys() and "LN" in chromosome.keys() and \
                chromosome["SN"].startswith(chromosome_startswith):
            chromosome_coverage_pl[chromosome["SN"]] = np.zeros(chromosome["LN"], dtype=np.uint32)
            chromosome_coverage_mn[chromosome["SN"]] = np.zeros(chromosome["LN"], dtype=np.uint32)
            chrom_size_list.append((chromosome["SN"], chromosome["LN"]))
            chromosome_pre_accessible_dict[chromosome["SN"]] = []
            chromosome_reads_dict[chromosome["SN"]] = 0
            # chromosome_lambda_dict[chromosome["SN"]] = 0.0
            genome_size += chromosome["LN"]

    if "output_chrom_size" in kwargs.keys() and kwargs["output_chrom_size"]:
        csize_fn = os.path.join(output_dir, output_prefix + ".csize")
        with open(csize_fn, "w") as f:
            pd.DataFrame(chrom_size_list, columns=["Chromosome", "Size"]).to_csv(f,
                                                                                 sep="\t",
                                                                                 index=False,
                                                                                 header=False)

    for read in bam:
        if read.mapq < mapq_threshold or \
                (library_type == "pe_fr" and read.is_reverse) or \
                (library_type == "pe_rf" and not read.is_reverse):
            continue
        if not read.reference_name.startswith(chromosome_startswith):
            continue

        read_strand = "-" if read.is_reverse else "+"
        pos_5prime = read.reference_start if read_strand == "+" else read.reference_end - 1
        chromosome_pre_accessible_dict[read.reference_name].append((read.reference_start, read.reference_end))
        chromosome_reads_dict[read.reference_name] += 1
        if read_strand == "+" and read.reference_name.startswith(chromosome_startswith):
            chromosome_coverage_pl[read.reference_name][pos_5prime] += 1
        elif read_strand == "-" and read.reference_name.startswith(chromosome_startswith):
            chromosome_coverage_mn[read.reference_name][pos_5prime] += 1

    # convert data type if necessary
    result_pl = {}
    result_mn = {}
    for chrom in chromosome_coverage_pl.keys():
        pl_chrom_max = np.max(chromosome_coverage_pl[chrom])
        mn_chrom_max = np.max(chromosome_coverage_mn[chrom])
        for k, v in enumerate(NP_DT_RANGE):
            if v > pl_chrom_max:
                chromosome_coverage_pl[chrom] = chromosome_coverage_pl[chrom].astype(NP_DT_NAME[k], copy=False)
                fn = os.path.join(output_dir, output_prefix + "_pl_%s" % chrom)
                np.save(fn, chromosome_coverage_pl[chrom])
                result_pl[chrom] = fn+".npy"
                break
        for k, v in enumerate(NP_DT_RANGE):
            if v > mn_chrom_max:
                chromosome_coverage_mn[chrom] = chromosome_coverage_mn[chrom].astype(NP_DT_NAME[k], copy=False)
                fn = os.path.join(output_dir, output_prefix + "_mn_%s" % chrom)
                np.save(fn, chromosome_coverage_mn[chrom])
                result_mn[chrom] = fn + ".npy"
                break
    return result_pl, result_mn


if __name__ == "__main__":
    pl, mn = get_read_signal("/local/storage/ly349/toolshed/BioQueue/workspace/2/uploads/PROseq/NG2017.bam",
                             loc_prime="R_5", chromosome_startswith="chr20", output_dir=".",
                             output_prefix="new_io_engine_test")
    pl_r, mn_r = get_read_signal("/local/storage/ly349/toolshed/BioQueue/workspace/2/uploads/PROseq/NG2017.bam",
                             loc_prime="R_5", chromosome_startswith="chr20", output_dir=".",
                             output_prefix="new_io_engine_test_r", reverse_complement=True)
    a = np.load(pl["chr20"])
    b = np.load(pl_r["chr20"])
    print("abc")