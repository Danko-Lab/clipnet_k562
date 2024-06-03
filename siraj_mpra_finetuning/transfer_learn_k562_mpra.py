import logging
import os
import sys
from pathlib import Path

import pandas as pd

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf
from tensorflow.keras.callbacks import CSVLogger, LearningRateScheduler
from tqdm.keras import TqdmCallback

sys.path.append("../../clipnet/")
import clipnet
import custom_loss
import mpra_gen
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


# Create data loaders for training and validation
ref_fp = "../data/mpra/k562_mpra_snps_ft_ref.fa.gz"
alt_fp = "../data/mpra/k562_mpra_snps_ft_alt.fa.gz"
mpra_fp = "../data/mpra/k562_allelic_mpra_snps.csv.gz"

chroms = (
    pd.read_csv("clipnet_data_fold_assignments.csv")
    .set_index("chrom")
    .to_dict()["fold"]
)

train_chroms = [k for k, v in chroms.items() if v not in [fold, fold % 9 + 1, 0]]
val_chroms = [k for k, v in chroms.items() if v == fold % 9 + 1]
print(f"Training on {train_chroms} and validating on {val_chroms}")

train_args = [
    ref_fp,
    alt_fp,
    mpra_fp,
    train_chroms,
    rnn_v10.batch_size,
]
val_args = [
    ref_fp,
    alt_fp,
    mpra_fp,
    val_chroms,
    rnn_v10.batch_size,
]
train_gen = mpra_gen.MPRAGen(*train_args)
val_gen = mpra_gen.MPRAGen(*val_args)

# Specify GPU usage
nn = clipnet.CLIPNET(n_gpus=1, use_specific_gpu=gpu)

# Load the reference and alternative models
outdir = Path(f"../models/clipnet_k562_mpra/f{fold}/")

ref_model = tf.keras.models.load_model(
    f"../models/clipnet_k562/fold_{fold}.h5", compile=False
)
alt_model = tf.keras.models.load_model(
    f"../models/clipnet_k562/fold_{fold}.h5", compile=False
)
ref_model.compile(
    optimizer=rnn_v10.optimizer(**rnn_v10.opt_hyperparameters),
    loss=rnn_v10.loss,
    loss_weights={"shape": 1, "sum": 1 / 500},
    metrics=rnn_v10.metrics,
)
alt_model.compile(
    optimizer=rnn_v10.optimizer(**rnn_v10.opt_hyperparameters),
    loss=rnn_v10.loss,
    loss_weights={"shape": 1, "sum": 1 / 500},
    metrics=rnn_v10.metrics,
)

# Create a new model that outputs the log2 fold change
logfc = tf.math.log((ref_model.output[1]) / (alt_model.output[1]))
mpra_net = tf.keras.Model(inputs=[ref_model.input, alt_model.input], outputs=logfc)
mpra_net.compile(
    optimizer=rnn_v10.optimizer(**rnn_v10.opt_hyperparameters),
    loss="mse",
    metrics=["mae", custom_loss.PearsonCorrelation],
)

# Compile the model
mpra_net_filepath = str(outdir.joinpath("mpra_net.h5"))
checkpoint = tf.keras.callbacks.ModelCheckpoint(
    mpra_net_filepath, verbose=0, save_best_only=True
)
early_stopping = tf.keras.callbacks.EarlyStopping(verbose=1, patience=10)
tqdm_callback = TqdmCallback(verbose=1, bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}")
csv_logger = CSVLogger(
    filename=outdir.joinpath("mpra_net.log"), separator=",", append=True
)

# Fit the model
fit_model = mpra_net.fit(
    x=train_gen,
    validation_data=val_gen,
    epochs=rnn_v10.epochs,
    steps_per_epoch=steps_per_epoch,
    verbose=0,
    callbacks=[
        checkpoint,
        early_stopping,
        tqdm_callback,
        csv_logger,
        LearningRateScheduler(warmup_lr),
    ],
)
