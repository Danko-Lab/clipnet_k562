"""
This file contains a number of functions and the class CGen (CLIPNET Generator) that
assist in loading data while training CLIPNET models.
"""

import logging
import os
import random
import sys

import numpy as np
import pandas as pd

sys.path.append("../../clipnet/")
import utils

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf
from tensorflow.keras.utils import Sequence


def load_data(
    data_fp: str,
    folds: list,
    cores=8,
    reverse_complement=False,
    model_type="regression",
):
    """
    Load a single, unfolded dataset. Use reverse_complement=True to load dataset
    reverse complemented.
    """
    # load data and check dimensions
    print(f"Loading data from {data_fp}.")
    data = pd.read_csv(data_fp)
    data = data[data.active_K562 == 1]
    # Filter data to only include the specified folds
    data = data[data["fold"].isin(folds)]
    # print(include[0])
    X = [
        utils.get_twohot_from_series(data["ref_seq"], cores=cores),
        utils.get_twohot_from_series(data["alt_seq"], cores=cores),
    ]
    if model_type == "regression":
        y = data.zscore_effect_size.values
    elif model_type == "classification":
        y = data.emVar_K562.values
    else:
        raise ValueError(
            f"Model type {model_type} must be 'regression' or 'classification'."
        )
    if reverse_complement:
        X = [utils.rc_twohot_het(x) for x in X]
    # output datasets
    return X, y


class MPRAGen(Sequence):
    def __init__(
        self,
        data_fp,
        folds,
        batch_size,
        in_window=1000,
        max_jitter=100,
        rc_augmentation=True,
        model_type="regression",
        swap_alleles=True,
        shuffle=True,
    ):
        self.batch_size = batch_size
        self.rc_augmentation = rc_augmentation
        self.swap_alleles = swap_alleles
        self.shuffle = shuffle
        self.model_type = model_type
        self.X, self.y = load_data(data_fp, folds, model_type=self.model_type)
        assert self.X[0].shape == self.X[1].shape, (
            "Reference and alternative sequence shapes must match."
            + f"{self.X[0].shape} != {self.X[1].shape}"
        )
        self.steps_per_epoch = len(self.y) // batch_size
        self.in_window = in_window
        self.trim = (self.X[0].shape[1] - in_window) // 2
        assert max_jitter <= self.trim, (
            "Max jitter must be less than or equal to the inferred max jitter."
            + f"{max_jitter} > {self.trim}"
        )
        self.max_jitter = max_jitter
        self.index = np.arange(len(self.y))
        self.on_epoch_end()

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
        if self.max_jitter > 0:
            jitter = random.randint(0, self.max_jitter)
            X = [x[:, jitter : self.in_window + jitter, :] for x in X]
        else:
            X = [x[:, self.trim : X.shape[1] - self.trim, :] for x in X]
        X = [tf.convert_to_tensor(x.copy(), dtype=tf.float32) for x in X]
        y = tf.convert_to_tensor(y.copy(), dtype=tf.float32)
        return X, y


def create_data_loader(generator):
    """
    Create a data loader for training and validation.
    """
    data_loader = tf.data.Dataset.from_generator(
        lambda: generator,
        output_signature=(
            (
                tf.TensorSpec(shape=(None, generator.in_window, 4), dtype=tf.float32),
                tf.TensorSpec(shape=(None, generator.in_window, 4), dtype=tf.float32),
            ),
            tf.TensorSpec(shape=(None, 1), dtype=tf.float32),
        ),
    )
    return data_loader, generator.steps_per_epoch
