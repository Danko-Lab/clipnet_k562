"""
This file contains a number of functions and the class CGen (CLIPNET Generator) that
assist in loading data while training CLIPNET models.
"""

import logging
import os
import random

import numpy as np
import pandas as pd
import pyfastx

import utils

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf


def load_data(
    ref_fp: str, alt_fp: str, mpra_fp: str, chroms=None, reverse_complement=False
):
    """
    Load a single, unfolded dataset. Use reverse_complement=True to load dataset
    reverse complemented.
    """
    # load data and check dimensions
    print(
        f"Loading sequence data from {ref_fp} and {alt_fp} and procap data from {mpra_fp}"
    )
    ref = utils.get_twohot_fasta_sequences(ref_fp)
    ref_chroms = [x.name.split(":")[0] for x in pyfastx.Fasta(ref_fp)]
    alt = utils.get_twohot_fasta_sequences(alt_fp)
    mpra = pd.read_csv(mpra_fp, sep="\t")
    y = np.log2(
        (mpra.mean_RNA_alt_K562 / mpra.mean_Plasmid_alt_K562)
        / (mpra.mean_RNA_ref_K562 / mpra.mean_Plasmid_ref_K562)
    )
    if chroms is not None:
        include = np.where(np.isin(ref_chroms, chroms))
        X = [ref[include], alt[include]]
        y = y[include]
    print("Successfully loaded data")
    # do rc_augmentation
    if reverse_complement:
        X = [utils.rc_twohot_het(x) for x in X]
    # output datasets
    return X, y


class MPRAGen(tf.keras.utils.Sequence):
    def __init__(
        self,
        ref_fp,
        alt_fp,
        mpra_fp,
        chroms,
        steps_per_epoch,
        batch_size,
        max_jitter=100,
        rc_augmentation=True,
        swap_alleles=True,
        shuffle=True,
    ):
        self.steps_per_epoch = steps_per_epoch
        self.batch_size = batch_size
        self.max_jitter = max_jitter
        self.rc_augmentation = rc_augmentation
        self.swap_alleles = swap_alleles
        self.shuffle = shuffle
        X, y = load_data(ref_fp, alt_fp, mpra_fp)
        self.X = X
        self.y = y
        self.index = np.arange(len(y))

    def __len__(self):
        """Denotes the number of batches per epoch"""
        return self.steps_per_epoch

    def on_epoch_end(self):
        """Shuffles fold indexes on start and after each epoch."""
        if self.shuffle:
            np.random.shuffle(self.index)

    def __getitem__(self, index):
        """
        Gets a batch of data.
        """
        batch_indices = self.index[
            index * self.batch_size : (index + 1) * self.batch_size
        ]
        X = [x[batch_indices, :, :] for x in self.X]
        y = self.y[batch_indices]
        if self.rc_augmentation and random.random() > 0.5:
            X = [x[:, ::-1, ::-1] for x in X]
        if self.swap_alleles and random.random() > 0.5:
            X = [X[1], X[0]]
            y = -y
        if self.max_jitter > 0:
            jitter = random.randint(0, self.max_jitter)
            X = [x[:, jitter : x.shape[1] + jitter, :] for x in X]
        return X, y
