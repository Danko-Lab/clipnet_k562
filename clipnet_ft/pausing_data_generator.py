import numpy as np
import tensorflow as tf
import utils


class DataGenerator(tf.keras.utils.Sequence):
    """
    Data generator for loading training/validation data efficiently.

    Args
    ----------
    X : np.array
        Input data. Shape (n_samples, n_features, n_channels).
    y : np.array
        Output data. Shape (n_samples, 1).
    batch_size : int
        Batch size.
    rc_augmentation : bool
        Whether to perform reverse complement augmentation.
    shuffle : bool
        Whether to shuffle data after each epoch.
    """

    def __init__(self, X, y, batch_size=64, rc_augmentation=True, shuffle=True):
        # Check that X and y have the same length
        if X.shape[0] != y.shape[0]:
            raise ValueError(f"lengths: X={X.shape[0]}, y={y.shape[0]}.")
        self.X = X
        self.y = y
        self.indexes = np.arange(self.X.shape[0])
        self.batch_size = batch_size
        self.rc_augmentation = rc_augmentation
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        """Denotes the number of batches per epoch"""
        return int(np.ceil(self.X.shape[0] / self.batch_size))

    def on_epoch_end(self):
        """Shuffles fold indexes on start and after each epoch."""
        if self.shuffle:
            np.random.shuffle(self.indexes)

    def __getitem__(self, index):
        """
        Gets data at a given index.
        """
        batch_indices = self.indexes[
            index * self.batch_size : (index + 1) * self.batch_size
        ]
        X = self.X[batch_indices, :, :]
        y = self.y[batch_indices]
        if self.rc_augmentation and np.random.rand() > 0.5:
            X = utils.rc_twohot_het(X)
        return X, y
