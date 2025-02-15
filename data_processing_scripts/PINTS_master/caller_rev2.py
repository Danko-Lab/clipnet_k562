#!/usr/bin/env python
# coding=utf-8
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2019-06-11
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
import argparse
import collections
import gzip
import logging
import os
import sys
import warnings
from multiprocessing import Pool

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections

warnings.filterwarnings("error")
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
DEFAULT_PREFIX = str(os.getpid())
logging.basicConfig(format='%(name)s - %(asctime)s - %(levelname)s: %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO,
                    handlers=[
                        logging.FileHandler(os.path.join(os.getcwd(), '%s.log' % DEFAULT_PREFIX)),
                        logging.StreamHandler()
                    ])
logger = logging.getLogger("PINTS - Caller")

try:
    import numpy as np
    import pandas as pd
    import pysam
    import requests
    import scipy
    from io_engine import get_coverage, get_coverage_bw, get_read_signal, log_assert
    from pybedtools import BedTool
    from scipy.signal import find_peaks, peak_widths
    from scipy.stats import binom_test, nbinom, poisson, probplot, uniform

    # from statsmodels.base.model import GenericLikelihoodModel
    from stats_engine import (
        get_outlier_threshold,
        pval_dist,
        zip_cdf,
        zip_em,
        zip_moment_estimators,
    )
    from statsmodels.stats.multitest import multipletests
except ImportError as e:
    missing_package = str(e).replace("No module named '", "").replace("'", "")
    logger.error("Please install %s first!" % missing_package)
    sys.exit(-1)

housekeeping_files = []
COMMON_HEADER = ('chromosome', 'start', 'end', 'name', 'padj', 'strand', 'reads',
                 'pval', 'mu_0', 'pi_0', 'mu_1', 'pi_1', 'summit')


def handle_exception(exc_type, exc_value, exc_traceback):
    """
    Handler for exception

    Parameters
    ----------
    exc_type :
    exc_value :
    exc_traceback :

    Returns
    -------

    Refs
    ----
    https://stackoverflow.com/questions/6234405/logging-uncaught-exceptions-in-python/16993115#16993115
    """
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


# redirect exception message to log
sys.excepthook = handle_exception


def check_version():
    """
    Check version

    Returns
    -------

    """
    try:
        r = requests.get('https://raw.githubusercontent.com/liyao001/PINTS/master/version.txt')
        latest_version = r.text.strip()
        latest_version_lst = list(map(int, latest_version.split(".")))
        version_file = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'version.txt')
        if not os.path.exists(version_file):
            logger.warning("Cannot locate the version file, please make sure the installation is complete.")
            return None
        fh = open(version_file, "r")
        local_ver = fh.read().strip()
        logger.info("PINTS version: %s" % local_ver)
        local_version = list(map(int, local_ver.split(".")))
        fh.close()
        for k, sub_ver in enumerate(latest_version_lst):
            if local_version[k] < sub_ver:
                logger.info("A later version of PINTS is available (%s)" % latest_version)
                break
    except Exception as e:
        logger.warning(e)


def run_command(cmd, repress_log=False):
    """
    Run command

    Parameters
    ----------
    cmd : str
        command
    repress_log : bool
        When it's set to False, if the command failed, the log will not be wrote to logger.

    Returns
    -------
    stdout : str
        Stdout output
    stderr : str
        Stderr output for the child process
    return_code : int
        Exit status of the child process
    """
    from subprocess import PIPE, Popen
    p = Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    stderr = stderr.decode("utf-8")
    stdout = stdout.decode("utf-8")
    if not repress_log:
        if p.returncode != 0:
            logger.error("Failed to run command %s" % cmd)
    return stdout, stderr, p.returncode


def runtime_check():
    """
    Runtime check, make sure all dependent tools are callable

    Parameters
    ----------
    not_found_rc : int
        return code from bash (or something like bash)
        when the called software / command cannot be found
    Returns
    -------

    """
    import shutil
    if sys.platform == "win32":
        logger.warning("No test had performed on Windows, so it might be buggy.")
    dependent_tools = ("bgzip", "tabix", "bedtools")
    for tool in dependent_tools:
        full_path = shutil.which(tool)
        if full_path is None:
            logger.error("Required tool %s is not callable" % tool)
            exit(1)


def merge_intervals(intervals, distance=0):
    """
    Merge intervals

    Parameters
    ----------
    intervals : tuple/list
        List / tuple of interval tuples
    distance : int
        Maximum distance between features allowed for features to be merged.
        Default is 0. That is, overlapping and/or book-ended features are merged.

    Returns
    -------
    merged_intervals : list
        Tuple of merged intervals

    Refs
    ----
        https://www.geeksforgeeks.org/merging-intervals/
    """
    log_assert(distance >= 0, "distance need to be >= 0", logger)
    s = sorted(intervals, key=lambda t: t[0])
    m = 0
    for t in s:
        if t[0] > s[m][1] + distance:
            m += 1
            s[m] = t[:2]
        else:
            # consider intervals
            # ((6, 8), (1, 9), (2, 4), (4, 7))
            # if we don't add an extra check
            # the final result will be (1, 8) instead of (1, 9)
            if s[m][1] <= t[1]:
                s[m] = (s[m][0], t[1])
    return s[:m + 1]


def sliding_window(chromosome_coverage, window_size, step_size):
    """
    Generate sliding windows

    Parameters
    ----------
    chromosome_coverage : array-like
        0-based per base coverage array for a certain chromosome
    window_size : int
        Window size for scanning
    step_size : int
        Step size for scanning

    Yields
    ------
    window : int
        Read counts in this window
    start : int
        0-based start coordinate of this window
    end : int
        0-based end coordinate of this window
    """
    if step_size < 1:
        logger.error("step_size must >= 1")
        raise ValueError("step_size must >= 1")
    if len(chromosome_coverage) < 1:
        logger.error("chromosome_coverage must >= 1")
        raise ValueError("chromosome_coverage must >= 1")

    total_bins = np.floor(chromosome_coverage.shape[0] / step_size - window_size / step_size + 1).astype(
        int)
    start = 0
    end = window_size
    for _ in range(total_bins):
        window = np.sum(chromosome_coverage[start:end])
        yield window, (start, end)
        start += step_size
        end = start + window_size


def _atom_ler(bed_handler, chromosome, query_start, query_end, small_window_threshold=4, peak_in_bg_threshold=1):
    """
    Atom operation for local env refinement

    Parameters
    ----------
    bed_handler : pysam.TabixFile
        pysam.TabixFile handler to the bed file which stores information about all peaks
    chromosome : str
        Name of the chromosome / contig
    query_start : int
        0-based start coordinate
    query_end : int
        0-based end coordinate
    small_window_threshold : int
        Candidate peaks with lengths shorter than this value will be skipped
    peak_in_bg_threshold : float
        Candidate peaks with density higher than this value will be removed from the local environment

    Returns
    -------
    se_coords : list of tuples
        List of peak coordinates (relative, tuple)
    re_coords : list of tuples
        List of peak coordinates (absolute, tuple)
    """
    se_coords = []
    re_coords = []
    try:
        query_start = query_start if query_start >= 0 else 0
        for sub_peak in bed_handler.fetch(chromosome, query_start, query_end, parser=pysam.asTuple()):
            sp_start = int(sub_peak[1])
            sp_start = sp_start if sp_start >= 0 else 0
            sp_end = int(sub_peak[2])
            if sp_start < query_start:
                sp_start = query_start
            if sp_end > query_end:
                sp_end = query_end
            if sp_end - sp_start < small_window_threshold:
                continue
            peak_dens = float(sub_peak[4])

            if peak_dens >= peak_in_bg_threshold:  # and peak_dens < peak_similarity_threshold:
                a = sp_start - query_start
                b = sp_end - query_start
                se_coords.append((a, b))
                re_coords.append((sp_start, sp_end))

    except ValueError as e:
        logger.error(str(e) + "\n%s:%d-%d" % (chromosome, query_start, query_end))
    return se_coords, re_coords


def remove_peaks_in_local_env(bed_handler, chromosome, query_start_left, query_end_left, query_start_right,
                              query_end_right, small_window_threshold, peak_in_bg_threshold, coverage_info,
                              fdr_target, cache):
    """
    IQR-based local environment refinement

    Parameters
    ----------
    bed_handler : pysam.TabixFile
        pysam.TabixFile handler to the bed file which stores information about all peaks
    chromosome : str
        Name of the chromosome / contig
    query_start_left : int
        0-based start coordinate for left-side local environment
    query_end_left : int
        0-based end coordinate for left-side local environment
    query_start_right : int
        0-based start coordinate for right-side local environment
    query_end_right : int
        0-based end coordinate for right-side local environment
    small_window_threshold : int
        Candidate peaks with lengths shorter than this value will be skipped
    peak_in_bg_threshold : float
        Candidate peaks with density higher than this value will be removed from the local environment
    coverage_info : np.array
        array which stores coverage info
    fdr_target : float
        fdr target

    Returns
    -------
    local_env : array-like
        Corrected local environment
    n_real_peaks_corrected : int
        Number of real peaks detected in local environment
    """
    ler_count = 0
    bg_mus = []
    local_env_left = coverage_info[query_start_left:query_end_left]
    local_env_right = coverage_info[query_start_right:query_end_right]
    se_l, re_l = _atom_ler(bed_handler=bed_handler, chromosome=chromosome, query_start=query_start_left,
                           query_end=query_end_left, small_window_threshold=small_window_threshold,
                           peak_in_bg_threshold=peak_in_bg_threshold)
    se_r, re_r = _atom_ler(bed_handler=bed_handler, chromosome=chromosome, query_start=query_start_right,
                           query_end=query_end_right, small_window_threshold=small_window_threshold,
                           peak_in_bg_threshold=peak_in_bg_threshold)
    coord_offset = len(local_env_left)
    new_local_env = np.concatenate((local_env_left, local_env_right), axis=None)

    uncertain_se_l = []
    uncertain_re_l = []
    uncertain_se_r = []
    uncertain_re_r = []

    for k, (b, s) in enumerate(re_l):
        cache_key = "%d-%d" % (b, s)
        if cache_key in cache:
            if cache[cache_key] == 1:
                new_local_env[se_l[k][0]:se_l[k][1]] = -1
                ler_count += 1
        else:
            uncertain_re_l.append(re_l[k])
            uncertain_se_l.append(se_l[k])
    for k, (b, s) in enumerate(re_r):
        cache_key = "%d-%d" % (b, s)
        if cache_key in cache:
            if cache[cache_key] == 1:
                new_local_env[se_r[k][0] + coord_offset:se_r[k][1] + coord_offset] = -1
                ler_count += 1
        else:
            uncertain_re_r.append(re_r[k])
            uncertain_se_r.append(se_r[k])

    se_coords = []
    se_coords.extend(uncertain_se_l)
    se_coords.extend([(b + coord_offset, s + coord_offset) for b, s in uncertain_se_r])

    re_coords = []
    re_coords.extend(uncertain_re_l)
    re_coords.extend(uncertain_re_r)

    n_candidate = len(uncertain_se_l) + len(uncertain_se_r)
    if n_candidate > 3:
        for b, s in uncertain_se_l:
            le = np.copy(local_env_left)
            le[b:s] = -1
            background_window = np.concatenate((le, local_env_right), axis=None)
            try:
                mu_, pi_, _, _ = zip_em(background_window[background_window >= 0])
            except Exception as e:
                logger.error(e)
            bg_mus.append(mu_)
        for b, s in uncertain_se_r:
            # se_coords.append((b + coord_offset, s + coord_offset))
            le = np.copy(local_env_right)
            le[b:s] = -1
            background_window = np.concatenate((le, local_env_left), axis=None)
            try:
                mu_, pi_, _, _ = zip_em(background_window[background_window >= 0])
            except Exception as e:
                logger.error(e)
            bg_mus.append(mu_)
        outlier_t, _ = get_outlier_threshold(bg_mus)
        for k, v in enumerate(bg_mus):
            cache_key = "%d-%d" % (re_coords[k][0], re_coords[k][1])
            if v < outlier_t:
                ler_count += 1
                new_local_env[se_coords[k][0]:se_coords[k][1]] = -1
                cache[cache_key] = 1
            else:
                cache[cache_key] = 0
        return new_local_env[new_local_env >= 0], ler_count
    elif n_candidate > 0:
        uncertain_re_l.extend(uncertain_re_r)
        # se_l.extend(se_r)
        for k, (b, s) in enumerate(uncertain_re_l):
            cache_key = "%d-%d" % (re_coords[k][0], re_coords[k][1])
            background_window = np.concatenate((coverage_info[b - 2500:b],
                                                coverage_info[s:s + 2500]),
                                               axis=None)
            mu_pk_em, pi_pk_em, _, _ = zip_em(coverage_info[b:s])
            mu_bg_em, pi_bg_em, _, _ = zip_em(background_window)

            p_val_formal = 1 - zip_cdf(mu_pk_em, pi_pk_em, mu_bg_em)
            if p_val_formal < fdr_target:
                new_local_env[se_coords[k][0]:se_coords[k][1]] = -1
                ler_count += 1
                cache[cache_key] = 1
            else:
                cache[cache_key] = 0

        return new_local_env[new_local_env >= 0], ler_count
    else:
        return new_local_env[new_local_env >= 0], ler_count


def check_window(coord_start, coord_end, mu_peak, pi_peak, chromosome_coverage, peak_in_bg_threshold,
                 mu_bkg_minimum, sp_bed_handler, chromosome_name, fdr_target, cache, small_window_threshold=5,
                 flanking=(10000, 5000, 1000)):
    """
    Calculate p-value for a peak

    Parameters
    ----------
    coord_start : int
        0-based start coordinate
    coord_end : int
        0-based end coordinate
    mu_peak : float
        mu_mle of the peak
    pi_peak : float
        pi_mle of the peak
    chromosome_coverage : array-like
        0-based per base coverage array for a certain chromosome
    peak_in_bg_threshold : float
        Candidate peaks with density higher than this value will be removed from the local environment
    mu_bkg_minimum : float

    sp_bed_handler :

    chromosome_name :

    small_window_threshold : int
        Candidate peaks with lengths shorter than this value will be skipped
    flanking : tuple
        Lengths of local environment that this function will check

    Returns
    -------
    p_value : float
        p_value for the peak
    window_value: int
        read counts in this window
    mu_0 : float
        mu for local env
    pi_0: float
        pi for local env
    """
    selected_window = chromosome_coverage[coord_start:coord_end]
    window_value = np.sum(selected_window)
    if coord_end - coord_start < small_window_threshold \
            or window_value == 0:
        return 1., window_value, 0, 0, (0, 0, 0)
    flanking = np.asarray(flanking, dtype=int) // 2
    mus = []
    pis = []
    ler_counts = []
    cache = dict()
    for k, f in enumerate(flanking):
        # cache = dict()
        qsl = coord_start - f
        qel = coord_start
        qsl = qsl if qsl >= 0 else 0
        qsr = coord_end
        qer = coord_end + f
        bg, x = remove_peaks_in_local_env(bed_handler=sp_bed_handler, chromosome=chromosome_name, query_start_left=qsl,
                                          query_end_left=qel, query_start_right=qsr, query_end_right=qer,
                                          small_window_threshold=small_window_threshold,
                                          peak_in_bg_threshold=peak_in_bg_threshold, coverage_info=chromosome_coverage,
                                          fdr_target=fdr_target, cache=cache)

        mu_, pi_, _, _ = zip_em(bg)
        mus.append(mu_)
        pis.append(pi_)
        ler_counts.append(x)

    mu_0 = np.mean(mus)  # mus[index]
    pi_0 = np.mean(pis)  # pis[index]
    if mu_bkg_minimum is not None and mu_0 < mu_bkg_minimum:
        mu_0 = mu_bkg_minimum
    # mu_1, pi_1, _, _ = zip_em(selected_window)
    # if pi_peak >
    pvalue = 1 - zip_cdf(mu_peak, pi_peak, mu_0)
    # pvalue = poisson.sf(mu_peak, mu_0)
    if pvalue == 0:
        pvalue = 10e-16

    return pvalue, window_value, mu_0, pi_0, ler_counts


def cut_peaks(window, peak_height, peak_threshold, peak_distance, peak_prominence, peak_width, peak_wlen,
              peak_rel_height, donor_tolerance, receptor_tolerance, ce_trigger, max_distance=20):
    """
    Cut peaks from the given window

    Parameters
    ----------
    window : array-like
        Per base read counts / coverage
    peak_height : number or ndarray or sequence, optional
        Required height of peaks. Either a number or None.
    peak_threshold : number or ndarray or sequence, optional
        Required threshold of peaks, the vertical distance to its neighbouring samples.
    peak_distance : number, optional
        Required minimal horizontal distance in samples between neighbouring peaks.
    peak_prominence : number or ndarray or sequence, optional
        Required prominence of peaks
    peak_width : number or ndarray or sequence, optional
        Required width of peaks in samples.
    peak_wlen : int, optional
        Used for calculation of the peaks prominences, thus it is only used if one of the
        arguments prominence or width is given.
    peak_rel_height : float, optional
        Used for calculation of the peaks width, thus it is only used if width is given.
    donor_tolerance : float
        From sub peak seeking for merging, the new density should be larger than dt*prev_d
    receptor_tolerance : float
        From sub peak seeking being merged, the new density should be larger than rt*prev_d
    ce_trigger : int
        Sub peak narrower than cet will trigger receptor tolerance check
    Returns
    -------
    merged_intervals : list
        List of tuples of merged intervals [(start_1, end_1), ... , (start_n, end_n)]
    """
    peaks, _ = find_peaks(window, height=peak_height, threshold=peak_threshold,
                          distance=peak_distance, prominence=peak_prominence,
                          width=peak_width, wlen=peak_wlen,
                          rel_height=peak_rel_height)
    widths, cor_heights, starts, ends = peak_widths(window, peaks, rel_height=peak_rel_height)
    intervals = []
    for k, start in enumerate(starts):
        intervals.append((int(start), int(ends[k])))
    mi = merge_intervals(intervals=intervals, distance=1)

    candidates = []
    for m in mi:
        events = 0
        for i in range(m[0], m[1] + 1):
            events += window[i]
        candidates.append((m[0], m[1], events, events / (m[1] - m[0])))
    fwd_search = []
    rev_search = []

    # forward search
    for k, c in enumerate(candidates):
        if k < len(candidates) - 1:
            new_total = c[2] + candidates[k + 1][2]
            new_density = new_total / (candidates[k + 1][1] - c[0])
            if new_density >= donor_tolerance * c[3]:
                distance_check = c[1] - c[0] < ce_trigger or candidates[k+1][0] - c[1] > max_distance
                if distance_check and new_density < (receptor_tolerance * candidates[k + 1][3]):
                    continue
                merged = (c[0], candidates[k + 1][1], new_total, new_density)
                fwd_search.append(merged)
            else:
                fwd_search.append(c)
        else:
            fwd_search.append(c)
    # reverse search
    for k in range(len(candidates) - 1, -1, -1):
        c = candidates[k]
        if k > 0:
            new_total = c[2] + candidates[k - 1][2]
            new_density = new_total / (c[1] - candidates[k - 1][0])
            if new_density >= donor_tolerance * c[3]:
                distance_check = c[1] - c[0] < ce_trigger or c[0] - candidates[k-1][1] > max_distance
                if distance_check and new_density < (receptor_tolerance * candidates[k - 1][3]):
                    continue
                merged = (candidates[k - 1][0], c[1], new_total, new_density)
                rev_search.append(merged)
            else:
                rev_search.append(c)
        else:
            rev_search.append(c)
    fwd_search.extend(rev_search)
    final = merge_intervals(fwd_search, distance=1)
    return final


def check_window_chromosome(rc_file, output_file, strand_sign, chromosome_name, window_size, step_size,
                            read_counts_threshold, peak_height, peak_threshold, peak_distance, peak_prominence,
                            peak_width, peak_wlen, peak_rel_height, fdr_target, ce_donor, ce_receptor, ce_trigger):
    """
    Evaluate windows on a chromosome

    Parameters
    ----------
    rc_file : str
        Path to numpy saved read coverage info
    output_file : str
        Path to store outputs
    strand_sign : str
        Strand of windows
    chromosome_name : str
        Name of this chromosome
    window_size : int
        Window size for binning and finding windows with read counts
    step_size : int
        Step size for binning and finding ...
    read_counts_threshold : int
        Windows / peaks with read counts fewer than this value will be skipped
    peak_height : number or ndarray or sequence, optional
        Required height of peaks. Either a number or None.
    peak_threshold : number or ndarray or sequence, optional
        Required threshold of peaks, the vertical distance to its neighbouring samples.
    peak_distance : number, optional
        Required minimal horizontal distance in samples between neighbouring peaks.
    peak_prominence : number or ndarray or sequence, optional
        Required prominence of peaks
    peak_width : number or ndarray or sequence, optional
        Required width of peaks in samples.
    peak_wlen : int, optional
        Used for calculation of the peaks prominences, thus it is only used if one of the
        arguments prominence or width is given.
    peak_rel_height : float, optional
        Used for calculation of the peaks width, thus it is only used if width is given.
    fdr_target : float
        fdr target
    ce_donor : float
        From sub peak seeking for merging, the new density should be larger than dt*prev_d
    ce_receptor : float
        From sub peak seeking being merged, the new density should be larger than rt*prev_d
    ce_trigger : int
        Sub peak narrower than cet will trigger receptor tolerance check

    Returns
    -------
    result_df : pd.DataFrame
        Window bed in dataframe
    """
    global housekeeping_files
    small_window_threshold = 5
    # ler_cache = dict()
    per_base_cov = np.load(rc_file, allow_pickle=True)
    subpeak_bed = output_file.replace(".bed", "_subpeaks_%s.bed" % chromosome_name)
    bins = []
    for window, coord in sliding_window(per_base_cov, window_size=window_size, step_size=step_size):
        if window > read_counts_threshold:  # no reads in the bin
            bins.append((chromosome_name, coord[0], coord[1], window))

    logger.info("Before merging, there are %d windows on %s" % (len(bins), chromosome_name))
    tmp_df = pd.DataFrame(bins, columns=("chromosome", "start", "end", "reads"))
    tmp_df["name"] = "."
    tmp_df["strand"] = strand_sign
    tmp_df = tmp_df.loc[:, ("chromosome", "start", "end", "name", "reads", "strand")]
    if tmp_df.shape[0] == 0:  # no hit
        return None
    # merge windows in case peaks are split into different windows
    bed_obj = BedTool(tmp_df.to_csv(sep="\t", index=False, header=False), from_string=True)
    bed_obj = bed_obj.merge(c=(4, 5, 6), o=("distinct", "sum", "distinct"))
    merged_fn = output_file.replace(".bed", "_%s_merged_windows.bed" % chromosome_name)
    bed_obj.moveto(merged_fn)
    merged_windows = pd.read_csv(merged_fn,
                                 sep="\t", names=["Chromosome", "Start", "End", "Name", "Reads", "Strand"])
    logger.info("After merging, there are %d windows on %s" % (merged_windows.shape[0], chromosome_name))

    bins = []
    all_peak_mus = []
    spb_fh = open(subpeak_bed, "w")
    index = 1
    for nr, row in merged_windows.iterrows():
        sub_peaks = cut_peaks(per_base_cov[row["Start"]:row["End"]],
                              peak_height=peak_height,
                              peak_threshold=peak_threshold,
                              peak_distance=peak_distance,
                              peak_prominence=peak_prominence,
                              peak_width=peak_width,
                              peak_wlen=peak_wlen,
                              peak_rel_height=peak_rel_height,
                              donor_tolerance=ce_donor,
                              receptor_tolerance=ce_receptor,
                              ce_trigger=ce_trigger,
                              )
        for sp in sub_peaks:
            start = sp[0] + row["Start"]
            end = sp[1] + row["Start"]
            peak_region = per_base_cov[start:end]
            window_value = peak_region.sum()
            n_start_sizes = sum(peak_region > 0)  # if n start sizes is smaller than 3, then ZIP shouldn't be used
            peak_len = end - start
            if peak_len < small_window_threshold \
                    or window_value == 0 or n_start_sizes < 3:
                mu_peak = window_value / peak_len
                pi_peak = 0
            else:
                mu_peak, pi_peak, _, _ = zip_em(peak_region)
                all_peak_mus.append(mu_peak)
            summit_coord = start + np.argmax(peak_region)
            spb_fh.write("%s\t%d\t%d\t%s-%d\t%f\t%s\t%f\t%d\n" % (
                chromosome_name, start, end, chromosome_name, index, mu_peak, strand_sign, pi_peak, summit_coord))
            index += 1
    spb_fh.close()

    _, _, rc = run_command("bgzip -f %s; tabix %s" % (subpeak_bed, subpeak_bed + ".gz"))
    log_assert(rc == 0 and os.path.exists(subpeak_bed + ".gz"),
               "Failed to generate sub peak bed file (%s%s)" % (chromosome_name, strand_sign),
               logger)

    bed_handler = pysam.TabixFile(subpeak_bed + ".gz")
    if len(all_peak_mus) == 0:
        logger.warning("No non-trivial peak was detected on chromosome %s" % chromosome_name)
        peak_threshold = 1
        bkg_mu_threshold = 0
    else:
        peak_threshold = 1
        bkg_mu_threshold = np.quantile(all_peak_mus, 0.1)
        logger.info("Minimum mu in local environment %f" % bkg_mu_threshold)
        # logger.info("Threshold for peak candidate in local environment %f" % peak_threshold)
    global_cache = dict()
    with gzip.open(subpeak_bed + ".gz", "rt") as peak_obj:
        for peak in peak_obj:
            candidate_peak = peak.split("\t")
            peak_start = int(candidate_peak[1])
            peak_end = int(candidate_peak[2])
            peak_id = candidate_peak[3]
            peak_mu = float(candidate_peak[4])
            peak_pi = float(candidate_peak[6])
            peak_summit = int(candidate_peak[7])
            pval, wv, mu_bg, pi_bg, lerc = check_window(coord_start=peak_start, coord_end=peak_end, mu_peak=peak_mu,
                                                        pi_peak=peak_pi, chromosome_coverage=per_base_cov,
                                                        peak_in_bg_threshold=peak_threshold,
                                                        mu_bkg_minimum=bkg_mu_threshold, sp_bed_handler=bed_handler,
                                                        chromosome_name=chromosome_name,
                                                        fdr_target=fdr_target, cache=global_cache)  # , cache=ler_cache)
            if wv > read_counts_threshold:
                bins.append((chromosome_name, peak_start, peak_end, peak_id, pval, wv, mu_bg, pi_bg, peak_mu, peak_pi,
                             peak_summit, lerc[0], lerc[1], lerc[2]))

    result_df = pd.DataFrame(bins, columns=("chromosome", "start", "end", "name", "pval", "reads",
                                            "mu_0", "pi_0", "mu_1", "pi_1", "summit", "ler_1", "ler_2", "ler_3"))
    # result_df["name"] = "."
    result_df["strand"] = strand_sign
    result_df = result_df.loc[:, ("chromosome", "start", "end", "name", "pval", "strand", "reads",
                                  "mu_0", "pi_0", "mu_1", "pi_1", "summit", "ler_1", "ler_2", "ler_3")]
    return result_df


def independent_filtering(df, output_to=None, **kwargs):
    """
    Independent filtering

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe of peak candidates (bed-like format)
    output_to : str or None, optional
        Path for the output
    **kwargs :
        Keyword arguments, fdr_target (float, default 0.1),
                           adjust_method (str, default fdr_bh),
                           ind_filter_granularity (float, default 0.005)
                           and output_diagnostic_plot (bool, default True) is effective
    Returns
    -------
    final_df : pd.DataFrame
        DataFrame with adjusted p-values

    Refs
    ----
    Bourgon, R., Gentleman, R. & Huber, W.
         Independent filtering increases detection power for high-throughput experiments.
         PNAS 107, 9546–9551 (2010).
    """
    log_assert(isinstance(df, pd.DataFrame), "df needs to be pd.DataFrame", logger)
    log_assert("reads" in df.columns, "Please provide read counts (reads) in df", logger)
    quantiles_tested = []
    windows_remain = []
    windows_rejected = []
    adjusted_df = None
    select_probe = None
    read_counts_threshold = 0
    for quantile in np.arange(0, 1, kwargs["ind_filter_granularity"]):
        threshold = df.reads.quantile(quantile)
        tmp_probe = df["reads"] > threshold
        filtered_df = df.loc[tmp_probe, :].copy()
        try:
            mult_res = multipletests(filtered_df["pval"], alpha=kwargs["fdr_target"], method=kwargs["adjust_method"])
            read_out = mult_res[1]
        except ZeroDivisionError:
            # if we cannot run BH, then we use Bonferroni correction
            read_out = filtered_df["pval"] * filtered_df.shape[0]
        filtered_df["padj"] = read_out
        # filtered_df.drop(columns="pval", inplace=True)
        filtered_df["name"] = "."
        filtered_df = filtered_df.loc[:, ("chromosome", "start", "end", "name", "padj", "strand", "reads")]
        # ift_result.append((quantile, filtered_df.shape[0], sum(read_out < fdr_target)))
        quantiles_tested.append(quantile)
        windows_remain.append(filtered_df.shape[0])
        windows_rejected.append(sum(read_out < kwargs["fdr_target"]))

        logger.info("Independent filtering: percentile %f, rejection %d (%d)" % (quantile,
                                                                                 windows_rejected[-1],
                                                                                 df.shape[0]))

        if len(windows_rejected) > 0:
            if windows_rejected[-1] >= max(windows_rejected):
                adjusted_df = filtered_df
                select_probe = tmp_probe
                read_counts_threshold = threshold
        else:
            adjusted_df = filtered_df
            select_probe = tmp_probe
            read_counts_threshold = threshold
    final_df = df.copy()
    final_df["padj"] = final_df["pval"]
    final_df.loc[select_probe, "padj"] = adjusted_df["padj"]
    quantiles_tested_arr = np.asarray(quantiles_tested)
    windows_rejected_arr = np.asarray(windows_rejected)
    if kwargs["output_diagnostic_plot"] and output_to is not None:
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(5.5, 5.5))
        ax2 = ax1.twinx()
        ax1.scatter(quantiles_tested, windows_remain, facecolors="none", edgecolors="#0571b0",
                    label="Windows remain")
        ax2.scatter(quantiles_tested, windows_rejected, facecolors="none", edgecolors="#ca0020", label="Rejections")
        ax2.annotate("Read counts: %d" % read_counts_threshold,
                     xy=(quantiles_tested_arr[np.argmax(windows_rejected_arr)], np.max(windows_rejected_arr)))
        ax1.set_xlabel("Quantile of filter")
        ax1.set_ylabel("Windows remain", color="#0571b0")
        ax1.tick_params(axis='y', labelcolor="#0571b0")
        ax2.set_ylabel("Rejections", color="#ca0020")
        ax2.tick_params(axis='y', labelcolor="#ca0020")
        plt.tight_layout()
        plt.savefig(output_to, transparent=True, bbox_inches="tight")
        plt.close()
        logger.info("Diagnostic plot for independent filtering was wrote to %s" % output_to)
    return final_df


def peaks_single_strand(per_base_cov, output_file, window_size, step_size, strand_sign, **kwargs):
    """
    Calling peaks on one strand

    Parameters
    ----------
    per_base_cov : dict
        Per base cov for available chromosomes
    output_file : str
        Path of output files
    window_size : int
        Window size for binning and finding windows with read counts
    step_size : int
        Step size for binning and finding ...
    strand_sign : str
        Strand sign for the data
    **kwargs :

    Returns
    -------
    result_df : pd.DataFrame
        All peaks on the specific strand of all chromosomes
    """
    global housekeeping_files
    fn, ext = os.path.splitext(output_file)

    args = []
    for chrom, pbc_npy in per_base_cov.items():
        sub_peaks_name = output_file.replace(".bed", "_subpeaks_%s.bed" % chrom)
        merged_name = output_file.replace(".bed", "_%s_merged_windows.bed" % chrom)
        args.append((pbc_npy, output_file, strand_sign, chrom, window_size, step_size,
                     kwargs["read_counts_threshold"], kwargs["peak_height"],
                     kwargs["peak_threshold"], kwargs["peak_distance"], kwargs["peak_prominence"],
                     kwargs["peak_width"], kwargs["peak_wlen"], kwargs["peak_rel_height"], kwargs["fdr_target"],
                     kwargs["donor_tolerance"], kwargs["receptor_tolerance"], kwargs["ce_trigger"]))
        housekeeping_files.append(merged_name)
        housekeeping_files.append(sub_peaks_name + ".gz")
        housekeeping_files.append(sub_peaks_name + ".gz.tbi")

    if kwargs["thread_n"] == 1:
        # for debugging
        sub_dfs = []
        for arg_i in args:
            sub_dfs.append(
                check_window_chromosome(arg_i[0], arg_i[1], arg_i[2], arg_i[3], arg_i[4], arg_i[5], arg_i[6], arg_i[7],
                                        arg_i[8], arg_i[9], arg_i[10], arg_i[11], arg_i[12], arg_i[13], arg_i[14],
                                        arg_i[15], arg_i[16], arg_i[17])
            )
    else:
        with Pool(kwargs["thread_n"]) as pool:
            sub_dfs = pool.starmap(check_window_chromosome, args)

    sub_dfs = [sdf for sdf in sub_dfs if sdf is not None]
    log_assert(len(sub_dfs) > 0, "No signal found across all chromosomes!", logger)
    tmp_df = pd.concat(sub_dfs)

    if kwargs["output_diagnostic_plot"]:
        tmp_df.to_csv(output_file.replace(".bed", "_debug.csv"), index=False)
    big_peaks_probe = tmp_df.end - tmp_df.start > kwargs["small_peak_threshold"]
    small_peaks_probe = tmp_df.end - tmp_df.start <= kwargs["small_peak_threshold"]
    lamb_global = tmp_df.loc[big_peaks_probe, "mu_1"].quantile(0.75)
    logger.info("Lambda for small peaks: %f" % lamb_global)
    inflated_small_peaks = np.sum(small_peaks_probe)
    tmp_df["pval"] = tmp_df.apply(lambda x:
                                  x["pval"] if x["end"] - x["start"] > kwargs["small_peak_threshold"] else poisson.sf(
                                      x["reads"],
                                      lamb_global *
                                      (x["end"] - x["start"])),
                                  axis=1)
    corrected_small_peaks = np.sum(np.logical_and(tmp_df["pval"] < kwargs["fdr_target"],
                                                  small_peaks_probe))
    logger.info("Significant small peaks after correction: %d (%d)" % (corrected_small_peaks, inflated_small_peaks))
    if kwargs["output_diagnostic_plot"]:
        pval_dist(tmp_df.loc[tmp_df["end"] - tmp_df["start"] > kwargs["small_peak_threshold"], "pval"],
                  logger=logger,
                  output_diagnostic_plot=kwargs["output_diagnostic_plot"],
                  output_to=fn + "_broad_pval_hist.pdf")
        pval_dist(tmp_df.loc[tmp_df["end"] - tmp_df["start"] <= kwargs["small_peak_threshold"], "pval"],
                  logger=logger,
                  output_diagnostic_plot=kwargs["output_diagnostic_plot"],
                  output_to=fn + "_narrow_peaks_pval_hist.pdf")
        pval_dist(tmp_df["pval"],
                  logger=logger,
                  output_diagnostic_plot=kwargs["output_diagnostic_plot"],
                  output_to=fn + "_pval_hist.pdf")

    # stratified independent filtering
    tmp_df_sm = independent_filtering(tmp_df.loc[small_peaks_probe, :], output_to=fn + "_idpf_sm.pdf", **kwargs)
    tmp_df_bg = independent_filtering(tmp_df.loc[big_peaks_probe, :], output_to=fn + "_idpf_bg.pdf", **kwargs)
    result_df = pd.concat([tmp_df_sm, tmp_df_bg])

    result_df = result_df.loc[:, COMMON_HEADER]
    result_df.sort_values(by=['chromosome', 'start'], inplace=True)
    result_df.to_csv(output_file, sep="\t", index=False, header=False)
    return output_file


def merge_opposite_peaks(sig_peak_bed, peak_candidate_bed, divergent_output_bed, bidirectional_output_bed,
                         singleton_bed, fdr_target, **kwargs):
    """
    Merge peaks on the opposite strand and generate divergent peak pairs

    Parameters
    ----------
    sig_peak_bed : str
        Path to bed file which contains significant peaks
    peak_candidate_bed : str
        Path to bed file which contains all candidate peaks on the opposite strand
    divergent_output_bed : str
        Path to output which stores divergent peaks
    bidirectional_output_bed : str
        Path to output which stores bidirectional peaks (divergent / convergent)
    singleton_bed : str
        Path to output which stores significant peaks which failed to pair

    **kwargs :

    Returns
    -------

    """
    tbx = pysam.TabixFile(peak_candidate_bed)
    fh = open(sig_peak_bed, "r")
    div_fh = open(divergent_output_bed, "w")
    bid_fh = open(bidirectional_output_bed, "w")
    sfp_fh = open(singleton_bed, "w")  # singletons failed to pair
    for nr, line in enumerate(fh):
        items = line.strip().split("\t")
        start = int(items[1])
        end = int(items[2])
        current_summit = int(items[-1])
        # allow overlapping
        if items[5] == "+":
            query_start = start - kwargs["close_threshold"]
            query_start = query_start if query_start >= 0 else 0
            query_end = end
        else:
            query_start = start
            query_end = end + kwargs["close_threshold"]

        opposite_start = np.nan
        opposite_end = np.nan
        opposite_pval = np.nan
        opposite_qval = np.nan
        opposite_sum = 0
        opposite_starts = []
        opposite_ends = []
        opposite_qvals = []
        opposite_pvals = []
        opposite_vals = []
        opposite_summits = []
        # since windows on each strand have been merged,
        # so here I expect the following iter returns at
        # most two records
        try:
            query_start = query_start if query_start >= 0 else 0
            for hit in tbx.fetch(items[0], query_start, query_end, parser=pysam.asTuple()):
                hit_start = int(hit[1])
                hit_end = int(hit[2])
                hit_score = float(hit[4])
                hit_reads = float(hit[6])  # in case the read counts had been normed
                opposite_summit = int(hit[-1])
                if hit_start < query_start:
                    continue
                if hit_end - hit_start < 3:  # filter single peaks
                    continue
                opposite_starts.append(hit_start)
                opposite_ends.append(hit_end)
                opposite_qvals.append(hit_score)
                opposite_pvals.append(float(hit[7]))
                opposite_vals.append(hit_reads)
                opposite_summits.append(opposite_summit)
            if len(opposite_pvals) > 0:
                index = np.argmin(opposite_pvals)
                opposite_start = opposite_starts[index]
                opposite_end = opposite_ends[index]
                opposite_pval = opposite_pvals[index]
                opposite_qval = opposite_qvals[index]
                opposite_summit = int(opposite_summits[index])
                opposite_sum = sum(opposite_vals[:index + 1])
        except ValueError as err:
            logger.warning("Problematic region for %s:%d-%d\n%s" % (items[0],
                                                                    query_start,
                                                                    query_end, err))
        if opposite_start is np.nan:
            sfp_fh.write(line)
        else:
            items.extend((str(opposite_start), str(opposite_end), str(opposite_pval), str(opposite_sum)))
            coords = (int(items[1]), int(items[2]), opposite_start, opposite_end)

            if items[5] == "+":
                fwd_summit = current_summit
                rev_summit = opposite_summit
            else:
                fwd_summit = opposite_summit
                rev_summit = current_summit

            tre_start = min(coords)
            tre_end = max(coords)
            if opposite_qval < fdr_target:
                pairing_confidence = "Stringent(qval)"
            elif opposite_pval < fdr_target:
                pairing_confidence = "Stringent(pval)"
            else:
                pairing_confidence = "Relaxed"
            if tre_end - tre_start > kwargs["div_size_min"]:
                candidate_values = (items[0], str(tre_start), str(tre_end), ".", items[4], items[5],
                                    str(float(items[6]) + opposite_sum), items[1], items[2],
                                    str(opposite_start), str(opposite_end), pairing_confidence + "\n")
                bid_fh.write("\t".join(candidate_values))
                if fwd_summit - rev_summit >= kwargs["summit_dist_min"]:
                    div_fh.write("\t".join(candidate_values))
            else:
                sfp_fh.write(line)
    fh.close()
    bid_fh.close()
    div_fh.close()
    sfp_fh.close()


def housekeeping():
    """
    Delete intermediate files

    Returns
    -------

    """
    global housekeeping_files
    try:
        for hf in housekeeping_files:
            if os.path.exists(hf):
                os.remove(hf)
    except Exception as e:
        logger.warning(str(e) + " (%s)" % hf)


def show_parameter_info(input_bam, output_dir, output_prefix, thread_n, **kwargs):
    """
    Show parameters

    Parameters
    ----------
    input_bam : str
        Path to the input
    output_dir : str
        Path to the output dir
    output_prefix : str
        Output prefix
    thread_n : int
        Number of threads

    kwargs

    Returns
    -------

    """
    logger.info("Parameters")
    logger.info("input_bam: %s" % input_bam)
    logger.info("output_dir: %s" % output_dir)
    logger.info("output_prefix: %s" % output_prefix)
    logger.info("thread_n: %s" % thread_n)
    for k, v in kwargs.items():
        logger.info("%s: %s" % (k, v))


def peak_calling(input_bam, output_dir=".", output_prefix="pints", filters=[],
                 thread_n=1, **kwargs):
    """
    Peak calling wrapper

    Parameters
    ----------
    input_bam : str
        Path to the input bam file
    output_dir : str
        Path to write output
    output_prefix : str
        Prefix for all outputs
    filters : list or tuple
        List of keywords to filter chromosomes
    thread_n : int
        Max number of sub processes that can be created.
    kwargs :

    Returns
    -------

    """
    logger.info("Start")
    global housekeeping_files
    # safety check
    if input_bam == "bigwig":
        log_assert(os.path.exists(kwargs["bw_pl"]), "Cannot find bigwig file %s" % kwargs["bw_pl"], logger)
        log_assert(os.path.exists(kwargs["bw_mn"]), "Cannot find bigwig file %s" % kwargs["bw_mn"], logger)
    else:
        log_assert(kwargs["bam_parser"] is not None, "Please specify which type of experiment this data "
                                                     "was generated from with --exp-type", logger)
        log_assert(os.path.exists(input_bam), "Cannot find input bam file %s" % input_bam, logger)
    log_assert(os.path.exists(output_dir) and os.path.isdir(output_dir), "Cannot write to %s" % output_dir, logger)

    runtime_check()
    check_version()
    show_parameter_info(input_bam, output_dir, output_prefix, thread_n, **kwargs)

    prefix = os.path.join(output_dir, output_prefix)

    if input_bam != "bigwig":
        chromosome_coverage_pl, chromosome_coverage_mn = get_read_signal(input_bam=input_bam,
                                                                         loc_prime=kwargs["bam_parser"],
                                                                         reverse_complement=kwargs["seq_rc"],
                                                                         output_dir=output_dir,
                                                                         output_prefix=output_prefix,
                                                                         filters=filters,
                                                                         **kwargs
                                                                         )
    else:
        chromosome_coverage_pl, chromosome_coverage_mn = get_coverage_bw(output_dir=output_dir,
                                                                         output_prefix=output_prefix,
                                                                         **kwargs)

    for v in chromosome_coverage_pl.values():
        housekeeping_files.append(v)

    for v in chromosome_coverage_mn.values():
        housekeeping_files.append(v)

    # peak calling (IQR)
    pl_peaks_bed = peaks_single_strand(per_base_cov=chromosome_coverage_pl,
                                       output_file=prefix + "_pl.bed",
                                       strand_sign="+",
                                       thread_n=thread_n,
                                       **kwargs)

    mn_peaks_bed = peaks_single_strand(per_base_cov=chromosome_coverage_mn,
                                       output_file=prefix + "_mn.bed",
                                       strand_sign="-",
                                       thread_n=thread_n,
                                       **kwargs)

    run_command("bgzip -f %s; tabix %s" % (pl_peaks_bed, prefix + "_pl.bed.gz"))
    run_command("bgzip -f %s; tabix %s" % (mn_peaks_bed, prefix + "_mn.bed.gz"))

    pl_df = pd.read_csv(prefix + "_pl.bed.gz", sep="\t", header=None,
                        names=COMMON_HEADER)
    mn_df = pd.read_csv(prefix + "_mn.bed.gz", sep="\t", header=None,
                        names=COMMON_HEADER)
    # filter huge merged bins
    pl_df = pl_df.loc[pl_df["end"] - pl_df["start"] < kwargs["window_size_threshold"], :]
    mn_df = mn_df.loc[mn_df["end"] - mn_df["start"] < kwargs["window_size_threshold"], :]

    sig_pl_bins = pl_df.loc[pl_df["padj"] < kwargs["fdr_target"], :]
    # sig_pl_bins.loc[:, "name"] = "Peak" + sig_pl_bins.index.map(str)
    with open(prefix + "_sig_pl.bed", "w") as f:
        sig_pl_bins.to_csv(f, sep="\t", index=False, header=False)
    merge_opposite_peaks(prefix + "_sig_pl.bed", prefix + "_mn.bed.gz",
                         divergent_output_bed=prefix + "_sig_pl_divergent_peaks.bed",
                         bidirectional_output_bed=prefix + "_sig_pl_bidirectional_peaks.bed",
                         singleton_bed=prefix + "_sig_pl_singletons.bed",
                         **kwargs)

    sig_mn_bins = mn_df.loc[mn_df["padj"] < kwargs["fdr_target"], :]
    # sig_mn_bins["name"] = "Peak" + sig_mn_bins.index.map(str)
    with open(prefix + "_sig_mn.bed", "w") as f:
        sig_mn_bins.to_csv(f, sep="\t", index=False, header=False)
    merge_opposite_peaks(prefix + "_sig_mn.bed", prefix + "_pl.bed.gz",
                         divergent_output_bed=prefix + "_sig_mn_divergent_peaks.bed",
                         bidirectional_output_bed=prefix + "_sig_mn_bidirectional_peaks.bed",
                         singleton_bed=prefix + "_sig_mn_singletons.bed",
                         **kwargs)

    command = "cat %s %s | sort -k1,1 -k2,2n | bedtools merge -i - -c 12 -o distinct > %s" % (
        prefix + "_sig_pl_bidirectional_peaks.bed", prefix + "_sig_mn_bidirectional_peaks.bed",
        prefix + "_bidirectional_peaks.bed")
    run_command(command)
    command = "cat %s %s | sort -k1,1 -k2,2n | bedtools merge -i - -c 12 -o distinct> %s" % (
        prefix + "_sig_pl_divergent_peaks.bed", prefix + "_sig_mn_divergent_peaks.bed", prefix + "_divergent_peaks.bed")
    run_command(command)
    command = "cat %s %s | sort -k1,1 -k2,2n > %s" % (prefix + "_sig_pl_singletons.bed",
                                                      prefix + "_sig_mn_singletons.bed",
                                                      prefix + "_single_peaks.bed")
    run_command(command)

    # generate ENCODE narrow peak file
    cmd = "cat %s %s | sort -k1,1 -k2,2n | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$4,$11,$6,$7,$8,$5,$13-$2}' > %s" % (
        prefix + "_sig_pl.bed", prefix + "_sig_mn.bed", prefix + "_narrowPeak.bed")
    run_command(cmd)

    housekeeping_files.append(prefix + "_sig_pl.bed")
    housekeeping_files.append(prefix + "_sig_mn.bed")
    housekeeping_files.append(prefix + "_sig_pl_singletons.bed")
    housekeeping_files.append(prefix + "_sig_mn_singletons.bed")
    housekeeping_files.append(prefix + "_sig_pl_divergent_peaks.bed")
    housekeeping_files.append(prefix + "_sig_mn_divergent_peaks.bed")
    housekeeping_files.append(prefix + "_sig_pl_bidirectional_peaks.bed")
    housekeeping_files.append(prefix + "_sig_mn_bidirectional_peaks.bed")

    logger.info("Finished")
    logger.info("Divergent peaks were saved to %s" % prefix + "_divergent_peaks.bed")
    logger.info("Bidirectional peaks were saved to %s" % prefix + "_bidirectional_peaks.bed")
    logger.info("Significant peaks which failed to pair were saved to %s" % prefix + "_single_peaks.bed")
    logger.info("Logs were saved to %s" % DEFAULT_PREFIX + ".log")
    # delete intermediate files
    housekeeping()
    return (prefix + "_divergent_peaks.bed", prefix + "_bidirectional_peaks.bed", prefix + "_single_peaks.bed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Peak Identifier for Nascent Transcript Sequencing")
    group = parser.add_argument_group("IO")
    group.add_argument("bam_file", help="input bam file, if you want to use bigwig files, please leave this as bigwig")
    group.add_argument("save_to", help="save peaks to")
    group.add_argument("file_prefix", action="store", default=str(os.getpid()),
                       help="prefix to all intermediate files")
    group.add_argument("--bw-pl", action="store", dest="bw_pl",
                       type=str, required=False,
                       help="Bigwig for plus strand. If you want to use bigwig instead of BAM, "
                            "please set bam_file to bigwig")
    group.add_argument("--bw-mn", action="store", dest="bw_mn",
                       type=str, required=False,
                       help="Bigwig for minus strand. If you want to use bigwig instead of BAM, "
                            "please set bam_file to bigwig")
    group.add_argument("--exp-type", action="store", default="CoPRO", dest="bam_parser",
                       help="Type of experiment, acceptable values are: CoPRO/GROcap/GROseq/PROcap/PROseq, or if you "
                            "know the position of RNA ends which you're interested on the reads, you can specify "
                            "R_5, R_3, R1_5, R1_3, R2_5 or R2_3")
    group.add_argument("--reverse-complement", action="store_true", dest="seq_reverse_complement",
                       required=False, default=False,
                       help="Set this switch if reads in this library represent the reverse complement of nascent "
                            "RNAs, like PROseq")
    group.add_argument('-f', '--filters', action='store', type=str, nargs="*", default=[],
                       help="reads from chromosomes whose names contain any matches in filters will be ignored")

    group = parser.add_argument_group("Filtering")
    group.add_argument("--mapq-threshold", action="store", dest="mapq_threshold",
                       type=int, required=False, default=30, help="Minimum mapping quality")
    group.add_argument("--close-threshold", action="store", dest="close_threshold",
                       type=int, required=False, default=300,
                       help="Distance threshold for two peaks (on opposite strands) to be merged")
    group.add_argument("--window-size", action="store", dest="window_size",
                       type=int, required=False, default=100, help="size for sliding windows")
    group.add_argument("--max-window-size", action="store", dest="window_size_threshold",
                       type=int, required=False, default=2000, help="max size of divergent windows")
    group.add_argument("--step-size", action="store", dest="step_size",
                       type=int, required=False, default=100, help="step size for sliding windows")
    group.add_argument("--fdr-target", action="store", dest="fdr_target",
                       type=float, required=False, default=0.1, help="FDR target for multiple testing")
    group.add_argument("--read-counts-threshold", action="store", dest="read_counts_threshold",
                       type=int, required=False, default=0,
                       help="Threshold for a window to be considered as having read")
    group.add_argument("--adjust_method", action="store", dest="adjust_method",
                       type=str, required=False, default="fdr_bh", help="method for calculating adjusted p-vals")
    group.add_argument("--small-peak-threshold", action="store", dest="small_peak_threshold",
                       type=int, required=False, default=5,
                       help="Threshold for small peaks, peaks with width smaller than this value will be required "
                            "to run extra test")
    group.add_argument("--ind-filtering-granularity", action="store", dest="ind_filter_granularity",
                       type=int, required=False, default=0.005,
                       help="Granularity for independent filtering")

    group = parser.add_argument_group("Edge trimming")
    group.add_argument("--donor-tolerance", action="store", dest="donor_tolerance",
                       type=float, required=False, default=1.0, help="Donor tolerance in best score segments")
    group.add_argument("--receptor-tolerance", action="store", dest="receptor_tolerance",
                       type=float, required=False, default=0.1, help="Receptor tolerance in best score segments")
    group.add_argument("--ce-trigger", action="store", dest="ce_trigger",
                       type=int, required=False, default=3, help="Trigger for receptor tolerance checking")

    group = parser.add_argument_group("Peak properties")
    group.add_argument("--top-peak-threshold", action="store", dest="top_peak_threshold",
                       type=float, required=False, default=0.9,
                       help="Min size for a divergent peak")
    group.add_argument("--peak-height", action="store", dest="peak_height",
                       type=int, required=False, default=None,
                       help="the minimal peak height")
    group.add_argument("--peak-threshold", action="store", dest="peak_threshold",
                       type=int, required=False, default=None,
                       help="Required threshold of peaks, the vertical distance to its neighbouring samples.")
    group.add_argument("--peak-distance", action="store", dest="peak_distance",
                       type=int, required=False, default=None,
                       help="Required minimal horizontal distance (>= 1) in samples between neighbouring peaks.")
    group.add_argument("--peak-prominence", action="store", dest="peak_prominence",
                       type=int, required=False, default=None,
                       help="Required prominence of peaks.")
    group.add_argument("--peak-width", action="store", dest="peak_width",
                       type=int, required=False, default=None,
                       help="Required width of peaks in samples.")
    group.add_argument("--peak-wlen", action="store", dest="peak_wlen",
                       type=int, required=False, default=None,
                       help="Used for calculation of the peaks prominences, "
                            "thus it is only used if one of the arguments prominence or width is given.")
    group.add_argument("--peak-rel-height", action="store", dest="peak_rel_height",
                       type=int, required=False, default=1,
                       help="Used for calculation of the peaks width, thus it is only used if width is given.")
    group.add_argument("--div-size-min", action="store", dest="div_size_min",
                       type=int, required=False, default=0,
                       help="Min size for a divergent peak")
    group.add_argument("--summit-dist-min", action="store", dest="summit_dist_min",
                       type=int, required=False, default=0,
                       help="Min dist between two summit")

    group = parser.add_argument_group("Other")
    group.add_argument("--chromosome-start-with", action="store", dest="chromosome_startswith",
                       type=str, required=False, default="chr",
                       help="Only keep reads mapped to chromosomes with this prefix")
    group.add_argument("--dont-output-chrom-size", action="store_false", dest="output_chrom_size",
                       required=False, default=True,
                       help="Don't write chromosome dict to local folder (not recommended)")
    group.add_argument("--output-diagnostic-plot", action="store_true", dest="output_diagnostic_plot",
                       required=False, default=False,
                       help="Save diagnostic plots (independent filtering and pval dist) to local folder")
    group.add_argument("--thread", action="store", dest="thread_n",
                       type=int, required=False, default=1,
                       help="Number of max threads")
    args = parser.parse_args()

    peak_calling(args.bam_file, args.save_to, args.file_prefix,
                 close_threshold=args.close_threshold, window_size=args.window_size, step_size=args.step_size,
                 fdr_target=args.fdr_target, adjust_method=args.adjust_method, mapq_threshold=args.mapq_threshold,
                 read_counts_threshold=args.read_counts_threshold, small_peak_threshold=args.small_peak_threshold,
                 chromosome_startswith=args.chromosome_startswith, peak_height=args.peak_height,
                 peak_threshold=args.peak_threshold, peak_distance=args.peak_distance,
                 window_size_threshold=args.window_size_threshold, peak_prominence=args.peak_prominence,
                 peak_width=args.peak_width, peak_wlen=args.peak_wlen, peak_rel_height=args.peak_rel_height,
                 output_chrom_size=args.output_chrom_size, output_diagnostic_plot=args.output_diagnostic_plot,
                 ind_filter_granularity=args.ind_filter_granularity, thread_n=args.thread_n,
                 div_size_min=args.div_size_min, summit_dist_min=args.summit_dist_min,
                 top_peak_threshold=args.top_peak_threshold, bw_pl=args.bw_pl, bw_mn=args.bw_mn,
                 donor_tolerance=args.donor_tolerance, receptor_tolerance=args.receptor_tolerance,
                 ce_trigger=args.ce_trigger, bam_parser=args.bam_parser, seq_rc=args.seq_reverse_complement,
                 filters=args.filters)
