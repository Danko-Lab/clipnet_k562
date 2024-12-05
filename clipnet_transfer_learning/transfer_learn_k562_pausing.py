# python transfer_learn_k562_proseq.py $fold $gpu

import json
import logging
import math
import os
import sys
from pathlib import Path

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.callbacks import CSVLogger, LearningRateScheduler
from tqdm.keras import TqdmCallback

sys.path.append("../../clipnet/")
import cgen
import clipnet
import custom_loss
import rnn_v10

fold = int(sys.argv[1])
gpu = int(sys.argv[2])


def warmup_lr(epoch, lr):
    """
    Learning rate warmup schedule.
    """
    print(f"LEARNING RATE = {lr}")
    if epoch < 1:
        return lr / 10
    elif epoch == 1:
        return lr * 10
    else:
        return lr


outdir = Path(f"../models/clipnet_k562_pausing/f{fold}/")
with open(outdir.joinpath("dataset_params.json"), "r") as f:
    dataset_params = json.load(f)
steps_per_epoch = math.floor(
    sum(dataset_params["n_samples_per_train_fold"]) * 2 / rnn_v10.batch_size
)
steps_per_val_epoch = math.floor(
    sum(dataset_params["n_samples_per_val_fold"]) * 2 / rnn_v10.batch_size
)
train_args = [
    dataset_params["train_seq"],
    dataset_params["train_proseq"],
    steps_per_epoch,
    rnn_v10.batch_size,
    dataset_params["pad"],
]
val_args = [
    dataset_params["val_seq"],
    dataset_params["val_proseq"],
    steps_per_val_epoch,
    rnn_v10.batch_size,
    dataset_params["pad"],
]
train_gen = cgen.CGen(*train_args)
val_gen = cgen.CGen(*val_args)

# Load pretrained model
nn = clipnet.CLIPNET(n_gpus=1, use_specific_gpu=gpu)
pretrained_model = tf.keras.models.load_model(
    f"../models/clipnet_k562_proseq/fold_{fold}.h5", compile=False
)

# Create new model architecture
new_output = layers.BatchNormalization(
    name="new_batch_normalization",
)(
    layers.Dense(1, name="new_dense")(
        layers.GlobalAvgPool1D(name="new_global_avg_pool_1d")(
            pretrained_model.get_layer("max_pooling1d_2").output
        )
    )
)
new_model = tf.keras.models.Model(
    inputs=pretrained_model.input,
    outputs=new_output,
)

# Compile
new_model.compile(
    optimizer=rnn_v10.optimizer(**rnn_v10.opt_hyperparameters),
    loss="msle",
    metrics=custom_loss.corr,
)

# Create callbacks
model_filepath = str(outdir.joinpath("clipnet_k562_pausing.h5"))
cp = tf.keras.callbacks.ModelCheckpoint(model_filepath, verbose=0, save_best_only=True)
early_stopping = tf.keras.callbacks.EarlyStopping(verbose=1, patience=20)
training_time = clipnet.TimeHistory()
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
    epochs=rnn_v10.epochs,
    steps_per_epoch=steps_per_epoch,
    verbose=0,
    callbacks=[
        cp,
        early_stopping,
        training_time,
        tqdm_callback,
        csv_logger,
        LearningRateScheduler(warmup_lr),
    ],
)
