import json
import logging
import math
import os
import sys
from pathlib import Path

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf
from tensorflow.keras.callbacks import CSVLogger, LearningRateScheduler
from tqdm.keras import TqdmCallback

sys.path.append("/home2/ayh8/clipnet/")
import cgen
import clipnet
import rnn_v10

run = int(sys.argv[1])
fold = int(sys.argv[2])
gpu = int(sys.argv[3])


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


outdir = Path(
    f"/home2/ayh8/clipnet_k562/models/clipnet_subsample_k562_models/n5/run_{run}/f{fold}/"
)
outdir.mkdir(parents=True, exist_ok=True)
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
    dataset_params["train_procap"],
    steps_per_epoch,
    rnn_v10.batch_size,
    dataset_params["pad"],
]
val_args = [
    dataset_params["val_seq"],
    dataset_params["val_procap"],
    steps_per_val_epoch,
    rnn_v10.batch_size,
    dataset_params["pad"],
]
train_gen = cgen.CGen(*train_args)
val_gen = cgen.CGen(*val_args)
nn = clipnet.CLIPNET(n_gpus=1, use_specific_gpu=gpu)
fit_model = tf.keras.models.load_model(
    f"/home2/ayh8/clipnet_subsampling/models/n5_run{run}/fold_{fold}.h5", compile=False
)
fit_model.compile(
    optimizer=rnn_v10.optimizer(**rnn_v10.opt_hyperparameters),
    loss=rnn_v10.loss,
    loss_weights={"shape": 1, "sum": dataset_params["weight"]},
    metrics=rnn_v10.metrics,
)
model_filepath = str(outdir.joinpath(f"../fold_{fold}.h5"))
cp = tf.keras.callbacks.ModelCheckpoint(model_filepath, verbose=0, save_best_only=True)
early_stopping = tf.keras.callbacks.EarlyStopping(verbose=1, patience=20)
training_time = clipnet.TimeHistory()
tqdm_callback = TqdmCallback(verbose=1, bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}")
csv_logger = CSVLogger(
    filename=outdir.joinpath("transfer.log"),
    separator=",",
    append=True,
)
fit_model = fit_model.fit(
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