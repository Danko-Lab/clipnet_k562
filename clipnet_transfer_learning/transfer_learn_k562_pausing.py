# python transfer_learn_k562_proseq.py $fold $gpu

import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import utils
from learning_rate_schedules import warmup_lr

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import custom_loss
import tensorflow as tf
from pausing_data_generator import DataGenerator
from tensorflow.keras import layers
from tensorflow.keras.callbacks import CSVLogger, LearningRateScheduler
from tqdm.keras import TqdmCallback

fold = int(sys.argv[1])

# Set VRAM usage to growth:
gpus = tf.config.experimental.list_physical_devices("GPU")
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)

# Load pausing indices
pausing_index_files = [
    f"../../data/pausing_index/k562_pausing_index_G{i}.bed" for i in (1, 5, 6)
]
names = ("chrom", "start", "stop", "name", "y", "strand")
pausing_index = pd.concat(
    [pd.read_csv(fp, header=None, sep="\t", names=names) for fp in pausing_index_files]
)

# Load sequences
twohot = utils.get_twohot_fasta_sequences(
    "../../data/pausing_index/k562_pausing_index.fa.gz"
)
sequences = np.concatenate([twohot for _ in pausing_index_files])

# Partition data into training and validation sets
test_folds = [fold]
val_folds = [(fold + 1) % 9 + 1]
train_folds = [i for i in range(1, 10) if i not in test_folds + val_folds]
chrom_splits = pd.read_csv("../data_processing/data_fold_assignments.csv")
train_chroms = [
    row["chrom"] for i, row in chrom_splits.iterrows() if row["fold"] in train_folds
]
val_chroms = [
    row["chrom"] for i, row in chrom_splits.iterrows() if row["fold"] in val_folds
]

train_idx = pausing_index.chrom.isin(train_chroms)
val_idx = pausing_index.chrom.isin(val_chroms)
train_gen = DataGenerator(sequences[train_idx], pausing_index.y[train_idx].to_numpy())
val_gen = DataGenerator(sequences[val_idx], pausing_index.y[val_idx].to_numpy())

# Load pretrained model
pretrained_model = tf.keras.models.load_model(
    f"../models/clipnet_k562_proseq/fold_{fold}.h5", compile=False
)

# Create new model architecture
new_output = layers.Activation("relu", name="new_relu")(
    layers.BatchNormalization(
        name="new_batch_normalization",
    )(
        layers.Dense(1, name="new_dense")(
            layers.GlobalAvgPool1D(name="new_global_avg_pool_1d")(
                pretrained_model.get_layer("max_pooling1d_2").output
            )
        )
    )
)
new_model = tf.keras.models.Model(
    inputs=pretrained_model.input,
    outputs=new_output,
)

# Compile
optimizer = tf.keras.optimizers.Adam
opt_hyperparameters = {
    "learning_rate": 0.001,
    "beta_1": 0.9,
    "beta_2": 0.999,
    "epsilon": 1e-7,
}

new_model.compile(
    optimizer=optimizer(**opt_hyperparameters),
    loss="msle",
    metrics=custom_loss.tf_spearmanr,
)

# Create callbacks
outdir = Path(f"../models/clipnet_k562_pausing/f{fold}/")
model_filepath = str(outdir.joinpath("clipnet_k562_pausing.h5"))
cp = tf.keras.callbacks.ModelCheckpoint(model_filepath, verbose=0, save_best_only=True)
early_stopping = tf.keras.callbacks.EarlyStopping(verbose=1, patience=10)
tqdm_callback = TqdmCallback(verbose=1, bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}")
csv_logger = CSVLogger(
    filename=outdir.joinpath("clipnet_k562_pausing.log"),
    separator=",",
    append=True,
)

# Fit
fit_model = new_model.fit(
    x=train_gen,
    validation_data=val_gen,
    epochs=100,
    steps_per_epoch=len(train_gen),
    verbose=0,
    callbacks=[
        cp,
        early_stopping,
        tqdm_callback,
        csv_logger,
        LearningRateScheduler(warmup_lr),
    ],
)
