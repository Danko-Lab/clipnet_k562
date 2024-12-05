# python transfer_learn_k562_proseq.py $fold $gpu

import logging
import os
import sys
from pathlib import Path

import pandas as pd
import pyfastx
import tqdm
from learning_rate_schedules import warmup_lr

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.callbacks import CSVLogger, LearningRateScheduler
from tqdm.keras import TqdmCallback

sys.path.append("../../clipnet/")
import clipnet
import custom_loss
import rnn_v10

fold = int(sys.argv[1])
gpu = int(sys.argv[2])

# Load pausing indices
pausing_index_files = [
    f"../../data/pausing_index/k562_pausing_index_G{i}.bed" for i in (1, 5, 6)
]
names = ("chrom", "start", "stop", "name", "y", "strand")
pausing_index = pd.concat(
    [pd.read_csv(fp, header=None, sep="\t", names=names) for fp in pausing_index_files]
)

# Load sequences
ref_genome = pyfastx.Fasta("../../data/pausing_index/hg38.fa.gz")
subsequences = [
    ref_genome.fetch(
        row["chrom"], (row["start"] - 500 + 1, row["start"] + 500), strand=row["strand"]
    )
    for i, row in tqdm.tqdm(pausing_index.iterrows())
]


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
    metrics=custom_loss.corr_log,
)

# Create callbacks
outdir = Path(f"../models/clipnet_k562_pausing/f{fold}/")
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
